#!/usr/bin/env python
"""
  NAME:
    gridCERESLib.py   

  PURPOSE:
    Read, analyze, and plot data from the Level3 CERES data

  PYTHON VERSION:
    3.7.4

  MODULES:
    - Matplotlib
    - datetime
    - sys
    - numpy
    - netCDF4
    
"""

import numpy as np
import sys
import gzip
import importlib
from datetime import datetime
from dateutil.relativedelta import relativedelta
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
#import matplotlib.colors as color
import matplotlib.path as mpath
from matplotlib.colors import rgb2hex,Normalize
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from os import system
import glob
from scipy.stats import pearsonr,spearmanr
from sklearn.linear_model import HuberRegressor
from sklearn.preprocessing import StandardScaler
from scipy.signal import argrelextrema, find_peaks
# The commands module was discontinued for Python 3, so if the user
# is using python 2, import commands instead

sys.path.append('/home/bsorenson')
from python_lib import circle, plot_trend_line, nearest_gridpoint, \
    aerosol_event_dict, init_proj, plot_lat_circles

##!## Compute a circle in axes coordinates, which we can use as a boundary
##!## for the map. We can pan/zoom as much as we like - the boundary will be
##!## permanently circular.
##!#theta = np.linspace(0, 2*np.pi, 100)
##!#center, radius = [0.5, 0.5], 0.5
##!#verts = np.vstack([np.sin(theta), np.cos(theta)]).T
##!#circle = mpath.Path(verts * radius + center)

datacrs = ccrs.PlateCarree() 
#mapcrs = ccrs.NorthPolarStereo()
mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

var_dict = {
    'SWF': 'CERES_SW_TOA_flux___upwards',
    'LWF': 'CERES_LW_TOA_flux___upwards',
    'clr': 'Clear_layer_overlap_percent_coverages',
    'cld': 'Clear_layer_overlap_percent_coverages',
    'alb': 'CERES_broadband_surface_albedo'
}
  
max_dict = {
    'SWF': 250.,
    'LWF': 370.,
    'TOTAL': 560.,
    'clr': 100.,
    'cld': 100.,
    'alb': 1.
}

min_dict = {
    'SWF': 120.,
    'LWF': 300.,
    'TOTAL': 450.,
    'clr': 0.,
    'cld': 0.,
    'alb': 0.
}

##!## Find the gridpoint in the gridded lat/lon data that 
##!## corresponds to the station at slat and slon
##!## ---------------------------------------------------- 
##!#def nearest_gridpoint(slat, slon, grid_lat, grid_lon):
##!#    fun_c = np.maximum(np.abs(grid_lat - slat), \
##!#        np.abs(grid_lon - slon))
##!#    m_idx = np.where(fun_c == np.min(fun_c))
##!#    return m_idx
##!#
##!#def covariance(x,y):
##!#    avg_x = np.average(x)
##!#    avg_y = np.average(y)
##!#    N = len(x)
##!#    if(len(x)!=len(y)):
##!#        print("ERROR: Arrays are not the same size.\nArray x has len=",len(x),\
##!#              "\nArray y has len=",len(y))
##!#    else:
##!#        cov = (np.sum((x-avg_x)*(y-avg_y)))/(N-1)
##!#        return cov
##!#
##!#def correlation(x,y):
##!#    avg_x = np.average(x)
##!#    avg_y = np.average(y)
##!#    N = len(x)
##!#    if(len(x)!=len(y)):
##!#        print("ERROR: Arrays are not the same size.\nArray x has len=",len(x),\
##!#              "\nArray y has len=",len(y))
##!#    else:
##!#        cov = (np.sum((x-avg_x)*(y-avg_y)))/(N-1)
##!#        std_x = np.std(x)
##!#        std_y = np.std(y)
##!#        corr = cov/(std_x*std_y)
##!#        return(corr)

def ice_trend_calc(years,months,ice,avg_ice,index,str_month):
    interpx = years[np.where(months==index)]
    interper = np.poly1d(np.polyfit(interpx,ice,1)) 
    # Normalize trend by dividing by number of years
    trend = (interper(interpx[-1])-interper(interpx[0]))
    pcnt_change = (trend/avg_ice)*100.
    print(str_month+" trend: ",np.round(trend,3)," % month mean: ",np.round(pcnt_change,3))
    return trend

def trend_calc(years,months,ice,avg_ice,index,str_month):
    interpx = years[np.where(months==index)]
    # Find the number of decades being analyzed
    num_dec = 216./120.
    if(len(interpx)!=0):
        #interpx = years
        interper = np.poly1d(np.polyfit(interpx,ice,1)) 
        # Normalize trend by dividing by number of years
        trend = (interper(interpx[-1])-interper(interpx[0]))
        pcnt_change = (trend/avg_ice)*100.
        # Find the percent change per decade
        pcnt_chg_dec = pcnt_change/num_dec
        print(str_month+" trend:\t",np.round(trend,3),"\t% month mean: ",np.round(pcnt_change,3),"\t%/decade: ",np.round(pcnt_chg_dec,3))
        return trend
    else:
        return -99.

# 2021/01/28: modified function to allow for Aqua data
def readgridCERES(start_date,end_date,param,satellite = 'Aqua',minlat=60.5,\
                 calc_month = True,season = ''):
    global CERES_data 
    CERES_data = {}
  
    spring=False
    summer=False
    autumn=False
    winter=False
    sunlight=False
    if(season=='spring'):
        spring=True
    elif(season=='summer'):
        summer=True
    elif(season=='autumn'):
        autumn=True
    elif(season=='winter'):
        winter=True
    elif(season=='sunlight'):
        sunlight=True
   

    # Grab all the files
    if(satellite == 'Terra'):
        base_path = '/home/bsorenson/data/CERES/SSF_1Deg/monthly/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    else:
        base_path = '/home/bsorenson/data/CERES/SSF_1Deg/monthly/Aqua/CERES_SSF1deg-Month_Aqua-MODIS_Ed4.1_Subset_'
    #base_path = '/data/CERES/SSF_1Deg/monthly/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    #base_path = '/home/bsorenson/data/CERES/SSF_1Deg/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    total_list = sorted(glob.glob(base_path+'*.nc')) 

    # Loop over all files and find the ones that match with the desired times
    start_date = datetime.strptime(str(start_date),'%Y%m')
    end_date = datetime.strptime(str(end_date),'%Y%m')
    final_list = []
    for f in total_list:
        fdate = f.split('_')[-1][:6]
        fdate = datetime.strptime(str(fdate),'%Y%m')
        if((fdate >= start_date) & (fdate <= end_date)):
            if(spring==True):
                if((fdate.month >= 3) & (fdate.month < 6)):
                    final_list.append(f)
            elif(summer==True):
                if((fdate.month >= 6) & (fdate.month < 9)):
                    final_list.append(f)
            elif(autumn==True):
                if((fdate.month >= 9) & (fdate.month < 12)):
                    final_list.append(f)
            elif(winter==True):
                if((fdate.month == 12) | (fdate.month < 3)):
                    final_list.append(f)
            elif(sunlight==True):
                if((fdate.month > 3) & (fdate.month < 10)):
                    final_list.append(f)
            else:
                final_list.append(f)
    time_dim = len(final_list)

    data = Dataset(final_list[0],'r')
    lat_ranges = data['lat'][:].data
    lon_ranges = data['lon'][:].data
    #lat_ranges = np.arange(60.5,89.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)


    lat_indices = np.where(lat_ranges>=minlat)[0]
    print(lat_indices)

    # Grid the lat and data
    grid_lon, grid_lat = np.meshgrid(lon_ranges,lat_ranges[lat_indices])

    CERES_data['param'] = param
    CERES_data['data']   = np.zeros((time_dim,len(lat_indices),len(lon_ranges)))
    CERES_data['trends'] = np.zeros((len(lat_indices),len(lon_ranges)))
    CERES_data['dates'] = [] 
    CERES_data['parm_name'] = data.variables[param].standard_name 
    CERES_data['parm_unit'] = data.variables[param].units 
    CERES_data['lat'] = grid_lat
    CERES_data['lon'] = grid_lon
    CERES_data['month_fix'] = ''
    CERES_data['season']=season
    CERES_data['satellite'] = satellite
    data.close()

    # Loop over good files and insert data into dictionary
    i_count = 0
    for ff in final_list:
        print(ff)
        data = Dataset(ff,'r')
        CERES_data['data'][i_count,:,:] = data.variables[param][0,lat_indices,:]
        CERES_data['dates'].append(ff.split('_')[-1][:6])
    #    print(data.variables[param][0,lat_indices,:])
        data.close() 
        i_count+=1
     
    #start_date = datetime.strptime(str(start_date),'%Y%m')
    #end_date = datetime.strptime(str(end_date),'%Y%m') +timedelta(days=31)
    #
    #tempdate = data.variables['time'].units
    #orig_date = datetime.strptime(tempdate.split()[2]+' '+tempdate.split()[3],'%Y-%m-%d %H:%M:%S')

    if(calc_month == True):
        CERES_data = calcCERES_MonthClimo(CERES_data)

    return CERES_data

# For now, assume that the user only wants to grab a single file of daily data
# User could cram multiple months into a single file, if wanted
def readgridCERES_daily(start_date,end_date,param,satellite = 'Aqua',minlat=60.5,season='all'):
    CERES_data = {}
  
    lat_ranges = np.arange(minlat,90.5,1.0)
    lon_ranges = np.arange(0.5,360.5,1.0)

    # Grab all the files
    if(satellite == 'Terra'):
        base_path = '/home/bsorenson/data/CERES/SSF_1Deg/daily/Terra/CERES_SSF1deg-Day_Terra-MODIS_Ed4.1_Subset_'
    else:
        base_path = '/home/bsorenson/data/CERES/SSF_1Deg/daily/Aqua/CERES_SSF1deg-Day_Aqua-MODIS_Ed4.1_Subset_'
    #base_path = '/data/CERES/SSF_1Deg/monthly/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    #base_path = '/home/bsorenson/data/CERES/SSF_1Deg/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    total_list = sorted(glob.glob(base_path+'*.nc'))

    # Loop over all files and find the ones that match with the desired times
    start_date = datetime.strptime(str(start_date),'%Y%m%d')
    end_date = datetime.strptime(str(end_date),'%Y%m%d')
    final_list = []
    for f in total_list:
        fdate = f.split('_')[-1][:6]
        fdate = datetime.strptime(str(fdate),'%Y%m')
        if((fdate>=start_date) & (fdate<=end_date)):
            final_list.append(f)
    time_dim = len(final_list)

    lat_indices = np.where(lat_ranges>=minlat)[0]

    # Open up the file
    data = Dataset(final_list[0],'r')
    print(final_list[0])

    CERES_data['param'] = param
    CERES_data['data']   = np.zeros((data.variables['time'].size,len(lat_indices),len(lon_ranges)))
    CERES_data['trends'] = np.zeros((len(lat_indices),len(lon_ranges)))
    CERES_data['dates'] = [] 
    CERES_data['parm_name'] = data.variables[param].standard_name 
    CERES_data['parm_unit'] = data.variables[param].units 
    CERES_data['lat'] = lat_ranges[lat_indices]
    CERES_data['lon'] = lon_ranges
    CERES_data['month_fix'] = ''
    CERES_data['season']=season
    CERES_data['satellite'] = satellite

    split_date = data.variables['time'].units.split()[2].split('-')
    base_date = datetime(year = int(split_date[0]), month = int(split_date[1]),\
         day = int(split_date[2]))

    # Loop over good files and insert data into dictionary
    for ti in range(data.variables['time'].size):
        CERES_data['data'][ti,:,:] = data.variables[param][ti,lat_indices,:]
        new_date = base_date + relativedelta(days = data.variables['time'][ti])
        CERES_data['dates'].append(new_date.strftime("%Y%m%d"))
        ### print(data.variables[param][0,lat_indices,:])
    data.close() 
     
    #start_date = datetime.strptime(str(start_date),'%Y%m')
    #end_date = datetime.strptime(str(end_date),'%Y%m') +timedelta(days=31)
    #
    #tempdate = data.variables['time'].units
    #orig_date = datetime.strptime(tempdate.split()[2]+' '+tempdate.split()[3],'%Y-%m-%d %H:%M:%S')

    return CERES_data

# Data period is of format YYYYMMDDHH
def readgridCERES_hrly_grid(data_dt,param,satellite = 'Aqua',minlat=60.0,season='all'):
    print(data_dt)
    lat_ranges = np.arange(minlat,90.0,0.25)
    lon_ranges = np.arange(-180.0,180.0,0.25)

    # Grab all the files
    if(satellite == 'Terra'):
        base_path = '/home/bsorenson/data/CERES/SSF_Level2/Terra/'
    else:
        base_path = '/home/bsorenson/data/CERES/SSF_Level2/Aqua/modis_comp/'
    total_list = sorted(glob.glob(base_path+'CERES_SSF_*.nc'))

    # Convert the desired dt to a datetime object to use for finding the file
    # -----------------------------------------------------------------------
    day = False
    if(len(data_dt) == 10):
        str_fmt = "%Y%m%d%H"
        hour_adder = 1
        hour_subtracter = 1
    elif(len(data_dt) == 8):
        print("Daily average")
        day = True
        str_fmt = "%Y%m%d"
        hour_subtracter = 0
        hour_adder = 24
        day_adder = 'daily_'
    else:
        print("INVALID PLOT TIME")
        return
    
    dt_data_begin = datetime.strptime(data_dt,str_fmt) \
        - relativedelta(hours = hour_subtracter)
    dt_data_end   = datetime.strptime(data_dt,str_fmt) + \
         relativedelta(hours = hour_adder)

    print(dt_data_begin,dt_data_end)

    # Loop over the list of current CERES data files and find the one that  
    # corresponds to the desired datetime
    # --------------------------------------------------------------------
    good_list = []
    for tfile in total_list:
        split_file = tfile.strip().split('/')[-1].split('_')
        begin_date = datetime.strptime(split_file[-1][:10],"%Y%m%d%H")
        end_date   = datetime.strptime(split_file[-1][11:21],"%Y%m%d%H")
        if((dt_data_begin >= begin_date) & (dt_data_end <= end_date) | \
           (dt_data_begin == end_date) | \
           ((day == True) & (begin_date >= dt_data_begin) & \
           (end_date <= dt_data_end))):
            print("Found matching file",tfile)
            good_list.append(tfile)
            #work_file = tfile

    #total_list = [base_path+'CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc']
    #total_list = subprocess.check_output('ls '+base_path+'CERES_SSF_Aqua-XTRK_Edition4A_Subset_'+year+date+'*.nc',\
    #          shell=True).decode('utf-8').strip().split('\n')
    
    # Set up values for gridding the AI data
    lat_gridder = minlat * 4.


    base_date = datetime(year=1970,month=1,day=1)

    # Assume that only one file will match up
    # --------------------------------------- 
    #for fileI in range(len(good_list)):
    data = Dataset(good_list[0],'r')
    time  = data.variables['time'][:]
    lat   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
    lon   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
    sza   = data.variables['CERES_solar_zenith_at_surface'][:]
    vza   = data.variables['CERES_viewing_zenith_at_surface'][:]
    azm   = data.variables['CERES_relative_azimuth_at_surface'][:]
    lon[lon>179.99] = -360.+lon[lon>179.99]
    ##!#if(param == 'clr'):
    ##!#    # clear layer overlap indices:
    ##!#    #    0 - clear
    ##!#    #    1 - lower
    ##!#    #    2 - upper
    ##!#    #    3 - upper over lower
    ##!#    flux  = data.variables[var_dict[param]][:,0]
    ##!#elif(param == 'cld'):
    ##!#    flux  = data.variables[var_dict[param]][:,1] + \
    ##!#            data.variables[var_dict[param]][:,2] 
    ##!#else:
    ##!#    flux  = data.variables[var_dict[param]][:]
    lwf  = data.variables['CERES_LW_TOA_flux___upwards'][:]
    swf  = data.variables['CERES_SW_TOA_flux___upwards'][:]
    time  = data.variables['time'][:]
   
    data.close()
 
    # Loop over the values and rows
    print(swf.shape) 

    # Convert each pixel relative time to absolute
    total_times = np.array([base_date + relativedelta(days = ttime) \
        for ttime in time])

    plt.close('all') 
    fig1 = plt.figure()
    plt.plot(total_times, azm, label = 'azm')
    plt.plot(total_times, lat, label = 'lat')
    plt.plot(total_times, vza, label = 'vza')   
    plt.legend()
    plt.show()
   
    # Extract only the data in the time window 
    # ----------------------------------------
    test_time = total_times[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_swf  = swf[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lwf  = lwf[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lat  = lat[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lon  = lon[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_sza  = sza[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_vza  = vza[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_azm  = azm[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]

    plt.close('all') 
    fig1 = plt.figure()
    plt.plot(test_time, test_azm, label = 'azm')
    plt.plot(test_time, test_lat, label = 'lat')
    plt.plot(test_time, test_vza, label = 'vza')   
    plt.legend()
    plt.show()


    # Determine where the LAT peaks are located and separated by 181
    # --------------------------------------------------------------
    lat_neg_peaks, _ = find_peaks(-test_lat, distance = 120)
    lat_peak_diffs = lat_neg_peaks[1:] - lat_neg_peaks[:-1]


    # Extract the data within the full cycles
    # ---------------------------------------    
    keep_lat_peaks = lat_neg_peaks[np.where((lat_peak_diffs == 181) | \
        (lat_peak_diffs == 182))]
    print(lat_peak_diffs)
    if(len(keep_lat_peaks) == 0):
        print("ERROR: no data within the desired time window")
        return

    begin_idx = keep_lat_peaks[0]
    end_idx   = keep_lat_peaks[-1]

    test_time = test_time[begin_idx: end_idx]
    test_swf  = test_swf[begin_idx: end_idx]
    test_lwf  = test_lwf[begin_idx: end_idx]
    test_lat  = test_lat[begin_idx: end_idx]
    test_lon  = test_lon[begin_idx: end_idx]
    test_sza  = test_sza[begin_idx: end_idx]
    test_vza  = test_vza[begin_idx: end_idx]
    test_azm  = test_azm[begin_idx: end_idx]

    # Make grid arrays for the data
    # -----------------------------
    grid_swf  = np.full((300, 181), np.nan)
    grid_lwf  = np.full((300, 181), np.nan)
    grid_lat  = np.full((300, 181), np.nan)
    grid_lon  = np.full((300, 181), np.nan)
    grid_sza  = np.full((300, 181), np.nan)
    grid_vza  = np.full((300, 181), np.nan)
    grid_azm  = np.full((300, 181), np.nan)

    # Loop over the data and insert into the grids
    for ii in range(len(keep_lat_peaks) - 1):
        idx_diff = keep_lat_peaks[ii+1] - keep_lat_peaks[ii]
        print(idx_diff)
        if(idx_diff == 181):
            grid_swf[ii,:len(test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_lwf[ii,:len(test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_lat[ii,:len(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_lon[ii,:len(test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_sza[ii,:len(test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_vza[ii,:len(test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_azm[ii,:len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
        elif(idx_diff > 181):
            grid_swf[ii,:len(test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
            grid_lwf[ii,:len(test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
            grid_lat[ii,:len(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
            grid_lon[ii,:len(test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
            grid_sza[ii,:len(test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
            grid_vza[ii,:len(test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
            grid_azm[ii,:len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]

    # Remove any rows in the grid arrays with any nans
    # ------------------------------------------------
    checkers = np.array([np.count_nonzero(~np.isnan(grid_lon[ii,:])) for ii in range(grid_lon.shape[0])])
    print(checkers) 

    grid_swf = grid_swf[np.where(checkers == 181)[0],:]
    grid_lwf = grid_lwf[np.where(checkers == 181)[0],:]
    grid_lat = grid_lat[np.where(checkers == 181)[0],:]
    grid_lon = grid_lon[np.where(checkers == 181)[0],:]
    grid_sza = grid_sza[np.where(checkers == 181)[0],:]
    grid_vza = grid_vza[np.where(checkers == 181)[0],:]
    grid_azm = grid_azm[np.where(checkers == 181)[0],:]

    grid_swf  = np.ma.masked_invalid(grid_swf)
    grid_lwf  = np.ma.masked_invalid(grid_lwf)
    grid_swf  = np.ma.masked_where((grid_swf > 5000) | (grid_swf == 0), grid_swf)
    grid_lwf  = np.ma.masked_where((grid_lwf > 5000) | (grid_lwf == 0), grid_lwf)

    # Put the data into a structure
    # -----------------------------
    CERES_grid_hrly = {}
    CERES_grid_hrly['param']     = param
    #CERES_grid_hrly['parm_name'] = var_dict[param]
    CERES_grid_hrly['date']      = data_dt
    CERES_grid_hrly['satellite'] = satellite
    CERES_grid_hrly['swf']       = grid_swf
    CERES_grid_hrly['lwf']       = grid_lwf
    CERES_grid_hrly['lat']       = grid_lat
    CERES_grid_hrly['lon']       = grid_lon  
    CERES_grid_hrly['vza']       = grid_vza
    CERES_grid_hrly['sza']       = grid_sza  
    CERES_grid_hrly['azm']       = grid_azm

    return CERES_grid_hrly  


    ##!## Plot a map
    ##!#fig2 = plt.figure()
    ##!#mapcrs = ccrs.NorthPolarStereo()
    ##!#ax0 = fig2.add_subplot(1,1,1, projection = mapcrs)
    ##!#mesh = ax0.pcolormesh(grid_lon, grid_lat, np.ma.masked_invalid(grid_flux),transform = datacrs, shading = 'auto')
    ##!#cbar = plt.colorbar(mesh,ax = ax0, orientation='vertical',\
    ##!#    shrink = 0.8)
    ##!#ax0.coastlines(resolution = '50m')
    ##!#ax0.set_extent([-180,180,65,90], datacrs)
    ##!#plt.show()

# Data period is of format YYYYMMDDHH
def readgridCERES_hrly(data_dt,param,satellite = 'Aqua',minlat=60.0,season='all'):
    CERES_data_hrly = {}

    lat_ranges = np.arange(minlat,90.0,0.25)
    lon_ranges = np.arange(-180.0,180.0,0.25)

    # Grab all the files
    if(satellite == 'Terra'):
        base_path = '/home/bsorenson/data/CERES/SSF_Level2/Terra/'
    else:
        base_path = '/home/bsorenson/data/CERES/SSF_Level2/Aqua/'
    total_list = sorted(glob.glob(base_path+'CERES_SSF_*.nc'))

    # Convert the desired dt to a datetime object to use for finding the file
    # -----------------------------------------------------------------------
    day = False
    if(len(data_dt) == 10):
        str_fmt = "%Y%m%d%H"
        hour_adder = 1
        hour_subtracter = 1
    elif(len(data_dt) == 8):
        print("Daily average")
        day = True
        str_fmt = "%Y%m%d"
        hour_subtracter = 0
        hour_adder = 24
        day_adder = 'daily_'
    else:
        print("INVALID PLOT TIME")
        sys.exit()
    
    dt_data_begin = datetime.strptime(data_dt,str_fmt) 
    #    relativedelta(hours = hour_subtracter)
    dt_data_end   = datetime.strptime(data_dt,str_fmt) + \
         relativedelta(hours = hour_adder)

    print(dt_data_begin,dt_data_end)

    # Loop over the list of current CERES data files and find the one that  
    # corresponds to the desired datetime
    # --------------------------------------------------------------------
    good_list = []
    for tfile in total_list:
        split_file = tfile.strip().split('/')[-1].split('_')
        begin_date = datetime.strptime(split_file[-1][:10],"%Y%m%d%H")
        end_date   = datetime.strptime(split_file[-1][11:21],"%Y%m%d%H")
        if((dt_data_begin >= begin_date) & (dt_data_end <= end_date) | \
           (dt_data_begin == end_date) | \
           ((day == True) & (begin_date >= dt_data_begin) & \
           (end_date <= dt_data_end))):
            print("Found matching file",tfile)
            good_list.append(tfile)
            #work_file = tfile

    #total_list = [base_path+'CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc']
    #total_list = subprocess.check_output('ls '+base_path+'CERES_SSF_Aqua-XTRK_Edition4A_Subset_'+year+date+'*.nc',\
    #          shell=True).decode('utf-8').strip().split('\n')
    
    # Set up values for gridding the AI data
    lat_gridder = minlat * 4.


    swf_grid = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
    base_date = datetime(year=1970,month=1,day=1)
 
    for fileI in range(len(good_list)):
        # read in data directly from HDF5 files
        print(good_list[fileI])
        data = Dataset(good_list[fileI],'r')
        lat   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        lon   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        sza   = data.variables['CERES_solar_zenith_at_surface'][:]
        lon[lon>179.99] = -360.+lon[lon>179.99]
        if(param == 'clr'):
            # clear layer overlap indices:
            #    0 - clear
            #    1 - lower
            #    2 - upper
            #    3 - upper over lower
            flux  = data.variables[var_dict[param]][:,0]
        elif(param == 'cld'):
            flux  = data.variables[var_dict[param]][:,1] + \
                    data.variables[var_dict[param]][:,2] 
        else:
            flux  = data.variables[var_dict[param]][:]
        time  = data.variables['time'][:]
    
        # Loop over the values and rows
        print(flux.shape) 
        for i in range(flux.shape[0]):
            #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
            local_time = base_date + relativedelta(days = time[i])
            if((local_time >= dt_data_begin) & (local_time < dt_data_end)):
                if((flux[i] < 5000) and (flux[i] > 0) and (lat[i] >= minlat)):
                   #and (sza[i] < 80)):
                    # Skip over the really small flux values when doing
                    # daily averages
                    if(day and (sza[i] >= 77)):
                        continue 
                    index1 = int(np.floor(lat[i]*4 - lat_gridder))
                    index2 = int(np.floor(lon[i]*4 + 720.))
                    
                    if(index1 < 0): index1 = 0
                    if(index1 > 719): index1 = 719
                    if(index2 < 0): index2 = 0                                                                                            
                    if(index2 > 1439): index2 = 1439
             
                    try: 
                        swf_grid[index2, index1] = (swf_grid[index2,index1]*count[index2,index1] + flux[i])/(count[index2,index1]+1)
                    except IndexError:
                        print(lat[i],lon[i],index1,index2)
                    #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                    #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                    #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                    count[index2, index1] = count[index2,index1] + 1
        data.close()
    
    CERES_data_hrly['param'] = param
    CERES_data_hrly['parm_name'] = var_dict[param]
    CERES_data_hrly['data']   = swf_grid
    CERES_data_hrly['counts'] = count
    CERES_data_hrly['date'] = data_dt
    CERES_data_hrly['lat'] = lat_ranges
    CERES_data_hrly['lon'] = lon_ranges
    CERES_data_hrly['satellite'] = satellite
     
    #start_date = datetime.strptime(str(start_date),'%Y%m')
    #end_date = datetime.strptime(str(end_date),'%Y%m') +timedelta(days=31)
    #
    #tempdate = data.variables['time'].units
    #orig_date = datetime.strptime(tempdate.split()[2]+' '+tempdate.split()[3],'%Y-%m-%d %H:%M:%S')

    return CERES_data_hrly

# This function assumes the data is being read from the CERES_data dictionary 
def calcCERES_MonthClimo(CERES_data):

    # Set up arrays to hold monthly climatologies
    month_climo = np.zeros((12,CERES_data['data'].shape[1],CERES_data['data'].shape[2]))

    # Mask the monthly averages
    local_data = np.copy(CERES_data['data'][:,:,:])
    local_mask = np.ma.masked_where(local_data == -999., local_data)

    # Calculate monthly climatologies
    for m_i in range(12):
        month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::12,:,:],axis=0)
        print("Month: ",CERES_data['dates'][m_i][4:]) 
    ##!## Calculate monthly climatologies
    ##!## Assume all Aqua CERES data starting at 200207 is read in previously
    ##!#month_climo[0,:,:]  = np.nanmean(local_mask[6::12,:,:],axis=0)  # January 
    ##!#month_climo[1,:,:]  = np.nanmean(local_mask[7::12,:,:],axis=0)  # February
    ##!#month_climo[2,:,:]  = np.nanmean(local_mask[8::12,:,:],axis=0)  # March 
    ##!#month_climo[3,:,:]  = np.nanmean(local_mask[9::12,:,:],axis=0)  # April 
    ##!#month_climo[4,:,:]  = np.nanmean(local_mask[10::12,:,:],axis=0) # May 
    ##!#month_climo[5,:,:]  = np.nanmean(local_mask[11::12,:,:],axis=0) # June 
    ##!#month_climo[6,:,:]  = np.nanmean(local_mask[0::12,:,:],axis=0)  # July 
    ##!#month_climo[7,:,:]  = np.nanmean(local_mask[1::12,:,:],axis=0)  # August 
    ##!#month_climo[8,:,:]  = np.nanmean(local_mask[2::12,:,:],axis=0)  # September 
    ##!#month_climo[9,:,:]  = np.nanmean(local_mask[3::12,:,:],axis=0)  # October 
    ##!#month_climo[10,:,:] = np.nanmean(local_mask[4::12,:,:],axis=0)  # November
    ##!#month_climo[11,:,:] = np.nanmean(local_mask[5::12,:,:],axis=0)  # December

    # Insert data into dictionary
    CERES_data['MONTH_CLIMO'] = month_climo

    return CERES_data

def calc_avgCERES(CERES_dict,save=False,start_date='200101',end_date='201812',minlat=70.5,month_fix=False):
    
    print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
    print("\nAlbedo Data\n")

    if(month_fix==True):
        CERES_dict['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = CERES_dict['lat']
    lon_ranges = CERES_dict['lon']
    # Create array to hold monthly averages over region
    initial_avgs = np.full(len(CERES_dict['dates']),-9999.)
    initial_years  = np.zeros(len(initial_avgs))
    initial_months = np.zeros(len(initial_avgs))
    
    # Loop over times and calculate average at each time
    for time in range(len(CERES_dict['dates'])):
        #temp = np.copy(CERES_dict['data'][time,:,:])
        #temp[temp==-999] = np.nan
        initial_years[time] = int(CERES_dict['dates'][time][:4])
        initial_months[time] = int(CERES_dict['dates'][time][4:])

        good_indices = np.where(CERES_dict['data'][time,:,:]!=-999)
        if(len(good_indices[0])!=0):
            temp = np.copy(CERES_dict['data'][time,:,:])
            temp[temp==-999] = np.nan
            initial_avgs[time] = np.nanmean(temp)
        else:
            initial_avgs[time] = -420
        # Remove October and November if desired to match up with Pistore et al 2014
        if((month_fix==True) & ((initial_months[time]==10) | (initial_months[time]==11) | (initial_months[time]==2))):
            initial_avgs[time] = -420

    # Remove missing values
    good_indices = np.where(initial_avgs!=-420)
    years = initial_years[good_indices]
    months = initial_months[good_indices]
    regional_avgs = initial_avgs[good_indices]

    #tttt[tttt==-999]=np.nan
    #np.nanmean(tttt)
    jan_alb = initial_avgs[good_indices][np.where(months==1)]
    feb_alb = initial_avgs[good_indices][np.where(months==2)]
    mar_alb = initial_avgs[good_indices][np.where(months==3)]
    apr_alb = initial_avgs[good_indices][np.where(months==4)]
    may_alb = initial_avgs[good_indices][np.where(months==5)]
    jun_alb = initial_avgs[good_indices][np.where(months==6)]
    jul_alb = initial_avgs[good_indices][np.where(months==7)]
    aug_alb = initial_avgs[good_indices][np.where(months==8)]
    sep_alb = initial_avgs[good_indices][np.where(months==9)]
    oct_alb = initial_avgs[good_indices][np.where(months==10)]
    nov_alb = initial_avgs[good_indices][np.where(months==11)]
    dec_alb = initial_avgs[good_indices][np.where(months==12)]

    # Calculate the mean albedo for each month
    avg_jan = np.average(jan_alb)
    avg_feb = np.average(feb_alb)
    avg_mar = np.average(mar_alb)
    avg_apr = np.average(apr_alb)
    avg_may = np.average(may_alb)
    avg_jun = np.average(jun_alb)
    avg_jul = np.average(jul_alb)
    avg_aug = np.average(aug_alb)
    avg_sep = np.average(sep_alb)
    avg_oct = np.average(oct_alb)
    avg_nov = np.average(nov_alb)
    avg_dec = np.average(dec_alb)
  
    # Calculate the average albedo value over the period
    avg_alb = np.average(initial_avgs[good_indices]) 
    print("Average albedo (N of ",str(int(CERES_dict['lat'][0])),") over period: ",avg_alb) 

    # Deseasonalize the data
    deseasonal_alb = np.copy(initial_avgs)
    #deseasonal_alb = np.copy(regional_avgs)
    
    deseasonal_alb[np.where(initial_months==1)]  = deseasonal_alb[np.where(initial_months==1)]  - avg_jan
    deseasonal_alb[np.where(initial_months==2)]  = deseasonal_alb[np.where(initial_months==2)]  - avg_feb
    deseasonal_alb[np.where(initial_months==3)]  = deseasonal_alb[np.where(initial_months==3)]  - avg_mar
    deseasonal_alb[np.where(initial_months==4)]  = deseasonal_alb[np.where(initial_months==4)]  - avg_apr
    deseasonal_alb[np.where(initial_months==5)]  = deseasonal_alb[np.where(initial_months==5)]  - avg_may
    deseasonal_alb[np.where(initial_months==6)]  = deseasonal_alb[np.where(initial_months==6)]  - avg_jun
    deseasonal_alb[np.where(initial_months==7)]  = deseasonal_alb[np.where(initial_months==7)]  - avg_jul
    deseasonal_alb[np.where(initial_months==8)]  = deseasonal_alb[np.where(initial_months==8)]  - avg_aug
    deseasonal_alb[np.where(initial_months==9)]  = deseasonal_alb[np.where(initial_months==9)]  - avg_sep
    deseasonal_alb[np.where(initial_months==10)] = deseasonal_alb[np.where(initial_months==10)] - avg_oct
    deseasonal_alb[np.where(initial_months==11)] = deseasonal_alb[np.where(initial_months==11)] - avg_nov
    deseasonal_alb[np.where(initial_months==12)] = deseasonal_alb[np.where(initial_months==12)] - avg_dec
    avg_deseasonal = np.average(deseasonal_alb[good_indices]) 
    print("Average deseasonal albedo over period: ",avg_deseasonal) 
    
    # Calculate the trend of the total data
    # ignore missing values
    interpx = np.arange(len(initial_years))
    total_interper = np.poly1d(np.polyfit(interpx[good_indices],initial_avgs[good_indices],1)) 
    # Normalize trend by dividing by number of years
    total_trend = (total_interper(interpx[good_indices][-1])-total_interper(interpx[good_indices][0]))
    print("Total albedo trend (200101 - 201812): ",np.round(total_trend,3))
    pcnt_change = (total_trend/avg_alb)*100.
    print("     % of average: ",pcnt_change)

    ##slope,intercept,r_value,p_value,std_err = scipy.stats.linregress(deseasonal_ext[good_indices],deseasonal_alb[good_indices])
    ##x_range = np.arange(min(deseasonal_ext[good_indices]),max(deseasonal_ext[good_indices]),0.01)
    ##newy = x_range*slope+intercept
    ##trend = newy[-1]-newy[0]

    # Calculate the trend of the deseasonalized data
    de_interpx = np.arange(len(initial_years))
    total_de_interper = np.poly1d(np.polyfit(de_interpx[good_indices],deseasonal_alb[good_indices],1)) 
    # Normalize trend by dividing by number of years
    total_de_trend = (total_de_interper(de_interpx[good_indices][-1])-total_de_interper(de_interpx[good_indices][0]))

    # Calculate the trend again using scipy
    # Also finds r_value and p_value
    slope,intercept,r_value,p_value,std_err = stats.linregress(de_interpx[good_indices],deseasonal_alb[good_indices])
    newy = de_interpx[good_indices]*slope+intercept
    test_total_de_trend = newy[-1]-newy[0] 

    print("Total deseasonal albedo trend (200101 - 201812): ",np.round(total_de_trend,5))
    print("            r_value                            : ",np.round(r_value,5))
    print("            p_value                            : ",np.round(p_value,5))
    pcnt_change = (total_de_trend/avg_alb)*100.
    print("     % of average: ",pcnt_change)
    
    # Calculate trend for each month data
    print("\nMonthly trends")
    trends = []
    trends.append(trend_calc(years,months,jan_alb,avg_jan,1,'January'))
    trends.append(trend_calc(years,months,feb_alb,avg_feb,2,'February'))
    trends.append(trend_calc(years,months,mar_alb,avg_mar,3,'March'))
    trends.append(trend_calc(years,months,apr_alb,avg_apr,4,'April'))
    trends.append(trend_calc(years,months,may_alb,avg_may,5,'May'))
    trends.append(trend_calc(years,months,jun_alb,avg_jun,6,'June'))
    trends.append(trend_calc(years,months,jul_alb,avg_jul,7,'July'))
    trends.append(trend_calc(years,months,aug_alb,avg_aug,8,'August'))
    trends.append(trend_calc(years,months,sep_alb,avg_sep,9,'September'))
    trends.append(trend_calc(years,months,oct_alb,avg_oct,10,'October'))
    trends.append(trend_calc(years,months,nov_alb,avg_nov,11,'November'))
    trends.append(trend_calc(years,months,dec_alb,avg_dec,12,'December'))
    
    
    fig0 = plt.figure(figsize=(9,9))
    ax1 = fig0.add_subplot(211)
    ax1.plot(interpx[good_indices],initial_avgs[good_indices],label='Average Albedo')
    ax1.plot(interpx[good_indices],total_interper(interpx[good_indices]),'--',label='Trend ('+str(np.round(total_trend,3))+'/Study Period)')
    ax1.set_title('Arctic (North of '+str(int(CERES_dict['lat'][0]))+' N) Average Albedo')
    #ax1.set_xlabel('Months After January 2001')
    ax1.set_xticks(np.arange(len(initial_years)+1)[::24])
    ax1.set_xticklabels(np.arange(2001,2020)[::2])
    ax1.set_ylabel('Albedo')
    ax1.legend()
    
    ax2 = fig0.add_subplot(212)
    ax2.plot(de_interpx[good_indices],deseasonal_alb[good_indices],label='Deseasonalized Average Albedo')
    ax2.plot(de_interpx[good_indices],total_de_interper(de_interpx[good_indices]),'--',label='Trend ('+str(np.round(total_de_trend,3))+'/Study Period)')
    ax2.text(de_interpx[good_indices][2],-0.02,"y = "+str(np.round(slope,4))+"x")
    ax2.text(de_interpx[good_indices][2],-0.03,"p = "+str(np.round(p_value,4)))
    ax2.set_title('Arctic (North of '+str(int(CERES_dict['lat'][0]))+' N) Deseasonalized Average Albedo')
    #ax2.set_xlabel('Months After January 2001')
    ax2.set_xticks(np.arange(len(initial_years)+1)[::24])
    ax2.set_xticklabels(np.arange(2001,2020)[::2])
    ax2.set_ylabel('Albedo Anomaly')
    ax2.legend()
   
    # Removing the monthly time series 
    ###ax2 = fig0.add_subplot(212)
    #### Plot the change for each month
    ####fig1 = plt.figure()
    ###ax2.plot(years[np.where(months==1)], jan_alb,label='January')
    ###ax2.plot(years[np.where(months==2)], feb_alb,label='February')
    ###ax2.plot(years[np.where(months==3)], mar_alb,label='March')
    ###ax2.plot(years[np.where(months==4)], apr_alb,label='April')
    ###ax2.plot(years[np.where(months==5)], may_alb,label='May')
    ###ax2.plot(years[np.where(months==6)], jun_alb,label='June')
    ###ax2.plot(years[np.where(months==7)], jul_alb,'--',label='July')
    ###ax2.plot(years[np.where(months==8)], aug_alb,'--',label='August')
    ###ax2.plot(years[np.where(months==9)], sep_alb,'--',label='September')
    ###ax2.plot(years[np.where(months==10)],oct_alb,'--',label='October')
    ###ax2.plot(years[np.where(months==11)],nov_alb,'--',label='November')
    ###ax2.plot(years[np.where(months==12)],dec_alb,'--',label='December')
    ###ax2.set_ylabel('Albedo')
    ###ax2.set_xticks(np.arange(2001,2019)[::2])
    ###ax2.legend(loc='upper center',bbox_to_anchor=(0.5,-0.05), ncol=6)
    if(save==True):
        outname = 'albedo_changes_with_deseason'+CERES_dict['month_fix']+'_'+str(int(CERES_dict['lat'][0]))+'.png'
        fig0.savefig(outname,dpi=300)
        print("Saved image ",outname)
   
    trends = np.array(trends)
    good_trends = trends[trends!=-99.] 
    fig3,ax = plt.subplots()
    labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    ax.plot(np.arange(2,len(good_trends)+2),good_trends)
    ax.set_xticks(np.arange(1,13))
    ax.set_xlim(1,13)
    ax.set_xticklabels(labels)
    ax.set_title('Trend in Monthly Arctic Albedo\n2001 - 2018')
    ax.set_ylabel('Albedo Trend (/Study Period)')
    plt.grid()
    if(save==True):
        outname = 'monthly_albedo_trends'+CERES_dict['month_fix']+'_'+str(int(CERES_dict['lat'][0]))+'.png'
        fig3.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()
    return initial_avgs,trends,deseasonal_alb

def calc_CERES_trend(CERES_dict,save=False,start_date='200012',end_date='201812',minlat=45.5,month_fix=False,\
                     adjusted=False):
    
    print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
    print("\nAlbedo Data\n")

    if(month_fix==True):
        CERES_dict['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = CERES_dict['lat']
    lon_ranges = CERES_dict['lon']
    ## Create array to hold monthly averages over region
    #initial_avgs = np.full(len(CERES_dict['dates']),-9999.)
    #initial_years  = np.zeros(len(initial_avgs))
    #initial_months = np.zeros(len(initial_avgs))
   
    # Loop over the lat and lon dimensions
    for xi in range(len(CERES_dict['lat'])):
        max_trend = -999
        min_trend = 999
        for yj in range(len(CERES_dict['lon'])):
            # Calculate the trend at the current box
            interp_data = CERES_dict['data'][:,xi,yj]
            good_indices = np.where(interp_data!=-9.99e+02)
            if(len(good_indices[0])==0):
                total_trend = -9999.
            else:
                interpx = np.arange(len(interp_data))
                #print(CERES_dict['data'][:,xi,yj][good_indices])
                total_interper = np.poly1d(np.polyfit( \
                    interpx[good_indices],\
                    CERES_dict['data'][:,xi,yj][good_indices],1)) 
                # Normalize trend by dividing by number of years
                total_trend = (total_interper(interpx[-1])-\
                               total_interper(interpx[0]))
                if(total_trend>max_trend):
                    max_trend = total_trend
                if(total_trend<min_trend):
                    min_trend = total_trend
            CERES_dict['trends'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    #plt.pcolormesh(plot_good_data,cmap=plt.cm.bwr,vmin=-50,vmax=50)
    #plt.colorbar(label='Albedo Trend')

    if(CERES_dict['param']=='toa_sw_clr_mon'):
        min_val=-40.
        max_val = 40.
        title_adder = 'TOA Clear-Sky SWF'
    elif(CERES_dict['param']=='toa_sw_all_mon'):
        min_val=-25.
        max_val = 25.
        title_adder = 'TOA All-Sky SWF'
    elif(CERES_dict['param']=='toa_sw_cld_mon'):
        min_val=-25.
        max_val = 25.
        title_adder = 'TOA Cloudy-Sky SWF'
    elif(CERES_dict['param']=='toa_lw_clr_mon'):
        min_val=-15.
        max_val = 15.
        title_adder = 'TOA Clear-Sky LWF'
    elif(CERES_dict['param']=='toa_lw_all_mon'):
        min_val=-10.
        max_val = 10.
        title_adder = 'TOA All-Sky LWF'
    elif(CERES_dict['param']=='toa_lw_cld_mon'):
        min_val=-10.
        max_val = 10.
        title_adder = 'TOA Cloudy-Sky LWF'
    elif(CERES_dict['param']=='toa_net_clr_mon'):
        min_val=-60.
        max_val = 60.
        title_adder = 'TOA Clear-Sky Net Flux'
    elif(CERES_dict['param']=='toa_net_all_mon'):
        min_val=-25.
        max_val = 25.
        title_adder = 'TOA All-Sky Net Flux'
    elif(CERES_dict['param']=='toa_net_cld_mon'):
        min_val=-25.
        max_val = 25.
        title_adder = 'TOA Cloudy-Sky Net Flux'
    elif(CERES_dict['param']=='toa_alb_clr_mon'):
        min_val=-0.15
        max_val = 0.15
        title_adder = 'TOA Clear-Sky Albedo'
    elif(CERES_dict['param']=='toa_alb_all_mon'):
        min_val=-0.06
        max_val = 0.06
        title_adder = 'TOA All-Sky Albedo'
    elif(CERES_dict['param']=='toa_alb_cld_mon'):
        min_val=-0.06
        max_val = 0.06
        title_adder = 'TOA Cloudy-Sky Albedo'
        
    colormap = plt.cm.bwr
    #coolwarm = plt.cm.coolwarm
    #ax = plt.axes(projection=ccrs.NorthPolarStereo())
    if(adjusted==True):
        fig1 = plt.figure(figsize=(8,5))
        ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=45.))
    else:
        fig1 = plt.figure()
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180,180,45,90],ccrs.PlateCarree())
    ax.gridlines()
    #ax = plt.axes(projection=ccrs.Miller())
    file_season = '_'+CERES_dict['season']
    title_season = ' '+CERES_dict['season']
    if(CERES_dict['season']=='all'):
        file_season=''
        title_season=''
        
    mesh = plt.pcolormesh(CERES_dict['lon'],CERES_dict['lat'],CERES_dict['trends'],\
            transform=ccrs.PlateCarree(),vmin=min_val,vmax=max_val,cmap=colormap)
    ax.set_title('Terra CERES '+title_adder+title_season.title()+' Trend\n'+\
                 CERES_dict['dates'][0]+' - '+CERES_dict['dates'][-1])
    ender = '.png'
    if(adjusted==True):
        ax.add_feature(cartopy.feature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
        axins = inset_axes(ax,width="5%",height="100%",loc='lower left',\
                    bbox_to_anchor=(1.05,0.,1,1),bbox_transform=ax.transAxes,\
                    borderpad=0)
        cbar = plt.colorbar(mesh,cax=axins,cmap=colormap,label=CERES_dict['parm_unit'])
        ax.set_xlim(-5170748.535086173,5167222.438879491)
        ax.set_ylim(-3913488.8763307533,3943353.899053069)
        ender = '_adjusted'+ender
    else:
        cbar = plt.colorbar(mesh,cmap=colormap,label=CERES_dict['parm_unit'])
    #cbar.set_label(parm_unit)
    ax.coastlines()
    if(save==True):
        outname = 'ceres_grid_terra_trend_'+CERES_dict['param']+'_'+\
                  CERES_dict['dates'][0]+'_'+CERES_dict['dates'][-1]+\
                  file_season+ender
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show() 
   
    #if(save==True):
    #    outname = 'monthly_albedo_trends'+CERES_dict['month_fix']+'_'+str(int(CERES_dict['lat'][0]))+'.png'
    #    fig3.savefig(outname,dpi=300)
    #    print("Saved image",outname)
    #else: 
    #    plt.show()
    #return initial_avgs,trends,deseasonal_alb

def make_scatter(CERES_dict_alb_clr,CERES_dict_sw_clr,\
                 CERES_dict_lw_clr,CERES_dict_net_clr):
    dot_size=8
    # Generate scatter plots
    title_text = 'December 2000 - December 2018'
    fig1 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_sw_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky SWF Trend')
    outname = 'ceres_clralb_clrsw_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)
    
    fig2 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_lw_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky LWF Trend')
    outname = 'ceres_clralb_clrlw_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)
    
    fig3 = plt.figure()
    plt.scatter(CERES_dict_alb_clr['trends'].flatten(),\
        CERES_dict_net_clr['trends'].flatten(),s=dot_size)
    plt.title(title_text)
    plt.xlabel('Terra CERES Clear-Sky Albedo Trend')
    plt.ylabel('Terra CERES Clear-Sky Net Flux Trend')
    outname = 'ceres_clralb_clrnet_trend_scatter.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)

def icePlots(infile,monthfix=False):

    print("\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =")
    print("\nIce Data\n")
    
    # Read data into arrays
    with open(infile,'r') as f:
        flines = f.readlines()
        initial_years  = np.zeros(len(flines)-1)
        initial_months = np.zeros(len(flines)-1)
        initial_extent = np.zeros(len(flines)-1)
        initial_area   = np.zeros(len(flines)-1)
        for i, line in enumerate(flines[1:]):
            templine = line.strip().split(',')
            initial_years[i] = templine[0]
            initial_months[i] = templine[1]
            initial_extent[i] = templine[4]
            initial_area[i] = templine[5]

    ### Remove winter months out of solidarity with albedo
    ##good_indices = np.where((initial_months<=2.)&(initial_months>=10.))
    ##years = initial_years[good_indices]
    ##months = initial_months[good_indices]
    ##extent = initial_extent[good_indices]
    years = initial_years
    months = initial_months
    extent = initial_extent
    
    jan_ext = extent[np.where(months==1)]
    feb_ext = extent[np.where(months==2)]
    mar_ext = extent[np.where(months==3)]
    apr_ext = extent[np.where(months==4)]
    may_ext = extent[np.where(months==5)]
    jun_ext = extent[np.where(months==6)]
    jul_ext = extent[np.where(months==7)]
    aug_ext = extent[np.where(months==8)]
    sep_ext = extent[np.where(months==9)]
    oct_ext = extent[np.where(months==10)]
    nov_ext = extent[np.where(months==11)]
    dec_ext = extent[np.where(months==12)]

    # Calculate the mean extent for each month
    avg_jan = np.average(jan_ext)
    avg_feb = np.average(feb_ext)
    avg_mar = np.average(mar_ext)
    avg_apr = np.average(apr_ext)
    avg_may = np.average(may_ext)
    avg_jun = np.average(jun_ext)
    avg_jul = np.average(jul_ext)
    avg_aug = np.average(aug_ext)
    avg_sep = np.average(sep_ext)
    avg_oct = np.average(oct_ext)
    avg_nov = np.average(nov_ext)
    avg_dec = np.average(dec_ext)
  
###    feby = feb_ext-np.full(len(feb_ext),avg_feb) 
###    sepy = sep_ext-np.full(len(sep_ext),avg_sep) 
###    figt = plt.figure()
###    plt.plot(feby,label='Feb anomaly')
###    plt.plot(sepy,label='Sep anomaly')
######    plt.plot(feb_ext,label='Feb')
######    plt.plot(np.full(len(feb_ext),avg_feb),label='Feb avg')
######    plt.plot(sep_ext,label='Sep')
######    plt.plot(np.full(len(sep_ext),avg_sep),label='Sep avg')
###    plt.legend()
###    plt.grid()
###    plt.show() 
    
    # Calculate the average extent value over the period
    avg_ext = np.average(extent) 
    #avg_ext = np.average(extent) 
    print("Average extent over period: ",avg_ext) 
   
    # Deseasonalize the data
    deseasonal_ext = np.copy(extent)
    
    deseasonal_ext[np.where(months==1)]  = deseasonal_ext[np.where(months==1)]  - avg_jan
    deseasonal_ext[np.where(months==2)]  = deseasonal_ext[np.where(months==2)]  - avg_feb
    deseasonal_ext[np.where(months==3)]  = deseasonal_ext[np.where(months==3)]  - avg_mar
    deseasonal_ext[np.where(months==4)]  = deseasonal_ext[np.where(months==4)]  - avg_apr
    deseasonal_ext[np.where(months==5)]  = deseasonal_ext[np.where(months==5)]  - avg_may
    deseasonal_ext[np.where(months==6)]  = deseasonal_ext[np.where(months==6)]  - avg_jun
    deseasonal_ext[np.where(months==7)]  = deseasonal_ext[np.where(months==7)]  - avg_jul
    deseasonal_ext[np.where(months==8)]  = deseasonal_ext[np.where(months==8)]  - avg_aug
    deseasonal_ext[np.where(months==9)]  = deseasonal_ext[np.where(months==9)]  - avg_sep
    deseasonal_ext[np.where(months==10)] = deseasonal_ext[np.where(months==10)] - avg_oct
    deseasonal_ext[np.where(months==11)] = deseasonal_ext[np.where(months==11)] - avg_nov
    deseasonal_ext[np.where(months==12)] = deseasonal_ext[np.where(months==12)] - avg_dec
    avg_deseasonal = np.average(deseasonal_ext) 
    print("Average deseasonal extent over period: ",avg_deseasonal) 

    # Calculate the trend of the total data
    interpx = np.arange(len(years))
    total_interper = np.poly1d(np.polyfit(interpx,extent,1)) 
    # Normalize trend by dividing by number of years
    total_trend = (total_interper(interpx[-1])-total_interper(interpx[0]))
    print("Total extent trend (200012 - 201812): ",np.round(total_trend,3))
    pcnt_change = (total_trend/avg_ext)*100.
    print("     % of average: ",pcnt_change)

    # Calculate the trend of the deseasonalized data
    interpx = np.arange(len(years))
    total_de_interper = np.poly1d(np.polyfit(interpx,deseasonal_ext,1)) 
    # Normalize trend by dividing by number of years
    total_de_trend = (total_de_interper(interpx[-1])-total_de_interper(interpx[0]))

    # Calculate the trend again using scipy
    # Also finds r_value and p_value
    slope,intercept,r_value,p_value,std_err = stats.linregress(interpx,deseasonal_ext)
    newy = interpx*slope+intercept
    test_total_de_trend = newy[-1]-newy[0] 

    print("Total deseasonal extent trend (200012 - 201812): ",np.round(total_de_trend,5))
    print("            r_value                            : ",np.round(r_value,5))
    print("            p_value                            : ",np.round(p_value,5))
    pcnt_change = (total_de_trend/avg_ext)*100.
    print("     % of average: ",pcnt_change)
    
    # Calculate trend for each month data
    print("\nMonthly trends")
    trends = []
    trends.append(trend_calc(years,months,jan_ext,avg_jan,1,'January'))
    trends.append(trend_calc(years,months,feb_ext,avg_feb,2,'February'))
    trends.append(trend_calc(years,months,mar_ext,avg_mar,3,'March'))
    trends.append(trend_calc(years,months,apr_ext,avg_apr,4,'April'))
    trends.append(trend_calc(years,months,may_ext,avg_may,5,'May'))
    trends.append(trend_calc(years,months,jun_ext,avg_jun,6,'June'))
    trends.append(trend_calc(years,months,jul_ext,avg_jul,7,'July'))
    trends.append(trend_calc(years,months,aug_ext,avg_aug,8,'August'))
    trends.append(trend_calc(years,months,sep_ext,avg_sep,9,'September'))
    trends.append(trend_calc(years,months,oct_ext,avg_oct,10,'October'))
    trends.append(trend_calc(years,months,nov_ext,avg_nov,11,'November'))
    trends.append(trend_calc(years,months,dec_ext,avg_dec,12,'December'))
   
    fig0 = plt.figure(figsize=(9,9))
    ax1 = fig0.add_subplot(211)
    ax1.plot(extent,label='Extent')
    ax1.plot(total_interper(interpx),'--',label='Trend ('+str(np.round(total_trend,3))+'/Study Period)')
    ax1.set_title('Arctic Sea Ice Extent')
    #ax1.set_xlabel('Months After January 2001')
    ax1.set_xticks(np.arange(len(years)+1)[::24])
    ax1.set_xticklabels(np.arange(2001,2020)[::2])
    ax1.set_ylabel('Extent (Mkm$^{2}$)')
    ax1.legend()

    ax2 = fig0.add_subplot(212)
    ax2.plot(deseasonal_ext,label='Deseasonalized Extent')
    ax2.plot(total_de_interper(interpx),'--',label='Trend ('+str(np.round(total_de_trend,3))+'/Study Period)')
    ax2.text(interpx[2],-1.00,"y = "+str(np.round(slope,4))+"x")
    ax2.text(interpx[2],-1.20,"p = "+str(np.round(p_value,4)))
    ax2.set_title('Arctic Deseasonalized Ice Extent')
    #ax2.set_xlabel('Months After January 2001')
    ax2.set_xticks(np.arange(len(years)+1)[::24])
    ax2.set_xticklabels(np.arange(2001,2020)[::2])
    ax2.set_ylabel('Extent Anomaly (Mkm$^{2}$)')
    ax2.legend()

    ###### Removing the monthly time series 
    ###ax2 = fig0.add_subplot(212)
    #### Plot the change for each month
    ####fig1 = plt.figure()
    ###ax2.plot(years[np.where(months==1)], jan_ext,label='January')
    ###ax2.plot(years[np.where(months==2)], feb_ext,label='February')
    ###ax2.plot(years[np.where(months==3)], mar_ext,label='March')
    ###ax2.plot(years[np.where(months==4)], apr_ext,label='April')
    ###ax2.plot(years[np.where(months==5)], may_ext,label='May')
    ###ax2.plot(years[np.where(months==6)], jun_ext,label='June')
    ###ax2.plot(years[np.where(months==7)], jul_ext,'--',label='July')
    ###ax2.plot(years[np.where(months==8)], aug_ext,'--',label='August')
    ###ax2.plot(years[np.where(months==9)], sep_ext,'--',label='September')
    ###ax2.plot(years[np.where(months==10)],oct_ext,'--',label='October')
    ###ax2.plot(years[np.where(months==11)],nov_ext,'--',label='November')
    ###ax2.plot(years[np.where(months==12)],dec_ext,'--',label='December')
    ###ax2.set_ylabel('Extent (Mkm$^{2}$)')
    ###ax2.set_xticks(np.arange(2001,2019)[::2])
    ###ax2.legend(loc='upper center',bbox_to_anchor=(0.5,-0.05), ncol=6)
    fig0.savefig('extent_changes_with_deseason.png',dpi=300)
    print("Saved image extent_changes_with_deseason.png") 
    ### Make a figure of the deseasonalized trends
    ##fig6 = plt.figure()
    ##plt.plot(deseasonal_ext)
    ##plt.ylabel('Extent Anomaly (Mkm$^{2}$)')
    ##plt.title('Deseasonalized Arctic Ice Extent')
    ##plt.show()
    
    
    fig3,ax = plt.subplots()
    labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    ax.plot(np.arange(1,13),trends)
    ax.set_xticks(np.arange(1,13))
    ax.set_ylim(-2.0,-0.5)
    ax.set_xticklabels(labels)
    ax.set_title('Trend in Monthly Sea Ice Extent\n2001 - 2018')
    ax.set_ylabel('Extent Trend (Mkm$^{2}$/Study Period)')
    plt.grid()
    fig3.savefig('monthly_extent_trends.png',dpi=300)
    print("Saved image monthly_extent_trends.png") 
    #plt.show()

    return extent,trends,months,years,deseasonal_ext

##!#def all_years(ice_avgs,alb_avgs,years,monthfix,minlat):
##!#    colors = np.arange(2001,2019) 
##!# 
##!#    cmap = cm.bwr
##!#    #norm = Normalize(vmin=vmin,vmax=vmax)
##!#    norm = matplotlib.colors.Normalize(vmin=2001,vmax=2018) 
##!#    mapper = ScalarMappable(norm=norm,cmap=cmap)
##!#    ccolors = ['red','darkorange','gold',\
##!#               'yellow','greenyellow','green',\
##!#               'springgreen','aquamarine','aqua',\
##!#               'deepskyblue','cornflowerblue','blue',\
##!#               'mediumslateblue','blueviolet','darkviolet',\
##!#               'purple','fuchsia','deeppink']
##!#    fig    = plt.figure(figsize=(8,7))
##!#    ax = fig.add_subplot(111)
##!#    #fig.set_figheight(8)
##!#    #fig.set_figwidth(10)
##!#    ax.set_xlim(3,17)
##!#    ax.set_ylim(0.3,0.65) 
##!#    patches = []
##!#    for xi in range(18):
##!#        pairs = []
##!#        for yi in range(12):
##!#            index = xi*12+yi
##!#            if(alb_avgs[index]!=-420):
##!#                pairs.append((ice_avgs[index],alb_avgs[index]))
##!#        color = mapper.to_rgba(2001+xi)
##!#        color2 = rgb2hex(color)
##!#        poly = Polygon(pairs,fill=False,edgecolor=color2) 
##!#        #poly = Polygon(pairs,fill=False,edgecolor=ccolors[xi]) 
##!#        plt.gca().add_patch(poly)
##!#        patches.append(poly) 
##!#    cax = fig.add_axes([0.91,0.11,0.025,0.77])
##!#    cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='vertical')
##!#    #p = PatchCollection(patches,alpha=0.4)
##!#    #p.set_array(np.array(colors))
##!#    #ax.add_collection(p)
##!#    #fig.colorbar(p,ax=ax)
##!#    ax.set_xlabel('Sea Ice Extent (Mkm$^{2}$)')
##!#    ax.set_ylabel('Clear-sky Albedo')
##!#    outname = 'yearly_scatter_ice_vs_alb'+CERES_dict['month_fix']+'_'+str(int(CERES_dict['lat'][0]))+'.png'
##!#    fig.savefig(outname,dpi=300)
##!#    print('Saved image ',outname)
##!#    #plt.show()

def plotCERES_hrly(pax, CERES_data_hrly, param, minlat=65, \
        vmin = None, vmax = None, title = '', label = '', \
        circle_bound = False, gridlines = True, grid_data = True, \
        zoom = True):

    if(vmin is None):
        vmin = min_dict[param.upper()]
    if(vmax is None):
        vmax = max_dict[param.upper()]

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.plasma
    ##!#if(pax is None):
    ##!#    if(minlat < 45):
    ##!#        mapcrs = ccrs.Miller()
    ##!#    else:
    ##!#        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    if(grid_data):
        plot_lon  = CERES_data_hrly['lon']
        plot_lat  = CERES_data_hrly['lat']
        if(param == 'total'):
            mask_flux = CERES_data_hrly['swf'] + CERES_data_hrly['lwf']
        else:
            mask_flux = CERES_data_hrly[param.lower()]
    else:
        plot_lat, plot_lon = np.meshgrid(CERES_data_hrly['lat'], \
            CERES_data_hrly['lon'])
        mask_flux = np.ma.masked_where(CERES_data_hrly['counts'] == 0, \
            local_data)
        local_data = np.copy(CERES_data_hrly[param])
 
    if(gridlines): 
        pax.gridlines()
    if(title == ''):
        title =  'CERES ' + param + ' ' + CERES_data_hrly['date']
    if(label == ''):
        label = param.upper() + ' [W/m$^{2}$]'
    pax.set_title(title)

    #plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
    mesh = pax.pcolormesh(plot_lon, plot_lat,mask_flux,transform = datacrs,\
        cmap = colormap, vmin = vmin, vmax = vmax, shading = 'auto')
    cbar = plt.colorbar(mesh,ax = pax, orientation='vertical',\
        extend = 'both')
    cbar.set_label(label,fontsize = 14, weight='bold')
    cbar.ax.tick_params(labelsize=14)

    pax.coastlines(resolution = '50m')
    #pax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)

    if(circle_bound):
        pax.set_boundary(circle, transform=pax.transAxes)

# Plot just a single swath on a 1-panel figure
#     skiprows - 
# NOTE: modified to work with readgridCERES_hrly_grid
# ---------------------------------------------------
def plotCERES_hrly_figure(date_str, param,  \
        only_sea_ice = False, minlat = 65., skiprows = None, \
        lat_circles = None, grid_data = True, zoom = False, \
        save = False):

    dt_date_str = datetime.strptime(date_str, "%Y%m%d%H%M")
    
    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(grid_data):
        CERES_data_hrly = readgridCERES_hrly_grid(date_str[:10], param, \
            satellite = 'Aqua', minlat = minlat)

        if(CERES_data_hrly is None):
            print("ERROR: no data returned from readgridCERES_hrly_grid")
            print("Quitting")
            return
    else:
        CERES_data_hrly = readgridCERES_hrly(date_str[:10], param, satellite = 'Aqua',\
            minlat = minlat, season = 'all')

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (6,6))
    #mapcrs = ccrs.NorthPolarStereo()
    mapcrs = init_proj(date_str)
    #mapcrs = ccrs.LambertConformal(central_longitude = -100.)
    ax0 = fig1.add_subplot(1,1,1, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    plotCERES_hrly(ax0, CERES_data_hrly, param, minlat = minlat, \
        vmin = None, vmax = None, title = '', label = '', \
        circle_bound = False, gridlines = False, grid_data = grid_data)

    #pax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    if(zoom):
        ax0.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)

    plt.suptitle(date_str)

    # ----------------------------------------------------
    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax0, lat_circles) 


    fig1.tight_layout()
    
    if(save == True):
        out_name = 'ceres_single_pass_' + CERES_data_hrly['param'].lower() + \
            '_' + CERES_data_hrly['date'] + '.png'
        plt.savefig(out_name,dpi=300)
        print('Saved image '+out_name)
    else:
        plt.show()


# This function automatically regenerates all known figures 
def plotCERES_AllPlots(CERES_dict,CERES_dict_summer,CERES_dict_winter,minlat=31., \
                       season='',trend_type='standard',sfc_type='all', \
                       start_date='200101',end_date='201512'):
    min_lat = minlat
    startdate = start_date
    enddate = end_date

    print("plotCERES whole_year all")
    plotCERES(CERES_dict,save=True,ptype='both',season='',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_Climo whole_year all")
    plotCERES_Climo(CERES_dict,save=True,ptype='both',season='',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_GlobalTrend whole_year all")
    plotCERES_GlobalTrend(CERES_dict,save=True,ptype='both',season='',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_LatTrend whole_year")
    plotCERES_LatTrend(CERES_dict,save=True,ptype='swf',season='',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat,check_lat=85.)
    print("plotCERES_LatClimo whole_year swf")
    plotCERES_LatClimo(CERES_dict,minlat=min_lat,ptype='swf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimo whole_year lwf")
    plotCERES_LatClimo(CERES_dict,minlat=min_lat,ptype='lwf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend whole_year swf")
    plotCERES_LatTrendTrend(CERES_dict,minlat=min_lat,ptype='swf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend whole_year lwf")
    plotCERES_LatTrendTrend(CERES_dict,minlat=min_lat,ptype='lwf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend whole_year swf")
    plotCERES_LatClimoTrend(CERES_dict,minlat=min_lat,ptype='swf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend whole_year lwf")
    plotCERES_LatClimoTrend(CERES_dict,minlat=min_lat,ptype='lwf',save=True,season='',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    if(CERES_dict['dtype']!='both'):
        print("plotCERES_ObCount whole_year swf")
        plotCERES_ObCount(CERES_dict,save=True,ptype='swf',season='',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_ObCount whole_year lwf")
        plotCERES_ObCount(CERES_dict,save=True,ptype='lwf',season='',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_LatObCount whole_year swf")
        plotCERES_LatObCount(CERES_dict,minlat=min_lat,ptype='swf',save=True,season='',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
        print("plotCERES_LatObCount whole_year lwf")
        plotCERES_LatObCount(CERES_dict,minlat=min_lat,ptype='lwf',save=True,season='',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)

    print("plotCERES summer all")
    plotCERES(CERES_dict_summer,save=True,ptype='both',season='summer',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_Climo summer all")
    plotCERES_Climo(CERES_dict_summer,save=True,ptype='both',season='summer',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_GlobalTrend summer all")
    plotCERES_GlobalTrend(CERES_dict_summer,save=True,ptype='both',season='summer',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_LatTrend summer")
    plotCERES_LatTrend(CERES_dict_summer,save=True,ptype='swf',season='summer',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat,check_lat=85.)
    print("plotCERES_LatClimo summer swf")
    plotCERES_LatClimo(CERES_dict_summer,minlat=min_lat,ptype='swf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimo summer lwf")
    plotCERES_LatClimo(CERES_dict_summer,minlat=min_lat,ptype='lwf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend summer swf")
    plotCERES_LatTrendTrend(CERES_dict_summer,minlat=min_lat,ptype='swf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend summer lwf")
    plotCERES_LatTrendTrend(CERES_dict_summer,minlat=min_lat,ptype='lwf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend summer swf")
    plotCERES_LatClimoTrend(CERES_dict_summer,minlat=min_lat,ptype='swf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend summer lwf")
    plotCERES_LatClimoTrend(CERES_dict_summer,minlat=min_lat,ptype='lwf',save=True,season='summer',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    if(CERES_dict_summer['dtype']!='both'):
        print("plotCERES_ObCount summer swf")
        plotCERES_ObCount(CERES_dict_summer,save=True,ptype='swf',season='summer',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_ObCount summer lwf")
        plotCERES_ObCount(CERES_dict_summer,save=True,ptype='lwf',season='summer',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_LatObCount summer swf")
        plotCERES_LatObCount(CERES_dict_summer,minlat=min_lat,ptype='swf',save=True,season='summer',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
        print("plotCERES_LatObCount summer lwf")
        plotCERES_LatObCount(CERES_dict_summer,minlat=min_lat,ptype='lwf',save=True,season='summer',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)


    print("plotCERES winter both")
    plotCERES(CERES_dict_winter,save=True,ptype='both',season='winter',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_Climo winter both")
    plotCERES_Climo(CERES_dict_winter,save=True,ptype='both',season='winter',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_GlobalTrend winter both")
    plotCERES_GlobalTrend(CERES_dict_winter,save=True,ptype='both',season='winter',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat)
    print("plotCERES_LatTrend winter")
    plotCERES_LatTrend(CERES_dict_winter,save=True,ptype='swf',season='winter',trend_type='standard',\
              start_date=startdate,end_date=enddate,minlat=min_lat,check_lat=85.)
    print("plotCERES_LatClimo winter swf")
    plotCERES_LatClimo(CERES_dict_winter,minlat=min_lat,ptype='swf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimo winter lwf")
    plotCERES_LatClimo(CERES_dict_winter,minlat=min_lat,ptype='lwf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend winter swf")
    plotCERES_LatTrendTrend(CERES_dict_winter,minlat=min_lat,ptype='swf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatTrendTrend winter lwf")
    plotCERES_LatTrendTrend(CERES_dict_winter,minlat=min_lat,ptype='lwf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend winter swf")
    plotCERES_LatClimoTrend(CERES_dict_winter,minlat=min_lat,ptype='swf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    print("plotCERES_LatClimoTrend winter lwf")
    plotCERES_LatClimoTrend(CERES_dict_winter,minlat=min_lat,ptype='lwf',save=True,season='winter',\
              trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
    if(CERES_dict_winter['dtype']!='both'):
        print("plotCERES_ObCount winter swf")
        plotCERES_ObCount(CERES_dict_winter,save=True,ptype='swf',season='winter',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_ObCount winter lwf")
        plotCERES_ObCount(CERES_dict_winter,save=True,ptype='lwf',season='winter',trend_type='standard',\
                  start_date=startdate,end_date=enddate,minlat=min_lat)
        print("plotCERES_LatObCount winter swf")
        plotCERES_LatObCount(CERES_dict_winter,minlat=min_lat,ptype='swf',save=True,season='winter',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)
        print("plotCERES_LatObCount winter lwf")
        plotCERES_LatObCount(CERES_dict_winter,minlat=min_lat,ptype='lwf',save=True,season='winter',\
                  trend_type='standard',sfc_type='all',start_date=startdate,end_date=enddate)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# NEW GRID DATA FUNCTIONS
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Calculation functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def calcCERES_grid_trend(CERES_data, month_idx, trend_type, minlat):

    if(month_idx == None):
        month_idx = 0
        index_jumper = 1
    else:
        month_adder = '_month'
        if(CERES_data['season'] == 'sunlight'):
            index_jumper = 6
        else:   
            index_jumper = 12

    lat_ranges = CERES_data['lat'][:,0]
    #lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    # Make copy of CERES_data array
    print(CERES_data['dates'][month_idx::index_jumper])
    local_data   = np.copy(CERES_data['data'][month_idx::index_jumper,:,:])
    #local_counts = np.copy(CERES_data['OB_COUNT'][month_idx::index_jumper,:,:])
    #local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where(((local_data == -999.) | \
        (CERES_data['lat'] < minlat)), local_data)
    ceres_trends = np.zeros(local_data.shape[1:])

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(lat_ranges)):
        for j in range(0,len(lon_ranges)):
            # Check the current max and min
            work_mask = local_mask[:,i,j][~local_mask[:,i,j].mask][0]
            if(len(work_mask) > 1):
                x_vals = np.arange(0,len(work_mask))
                # Find the slope of the line of best fit for the time series of
                # average data
                if(trend_type=='standard'): 
                    slope, intercept, r_value, p_value, std_err = \
                        stats.linregress(x_vals,work_mask)
                    ceres_trends[i,j] = slope * len(x_vals)
                else:
                    res = stats.theilslopes(work_mask, x_vals, 0.90)
                    ceres_trends[i,j] = res[0]*len(x_vals)
            else:
                print('no data')

    ceres_trends = np.ma.masked_where(((CERES_data['lat'] < minlat) | \
        (ceres_trends == -999.)), ceres_trends)

    return ceres_trends

# title is the plot title
# ptype is 'climo', 'trend'
def plotCERES_spatial(pax, plat, plon, pdata, ptype, ptitle = '', plabel = '', \
        vmin = None, vmax = None, colorbar_label_size = 16, minlat = 65.):

    if(vmin == None):
        vmin = np.nanmin(pdata)
    if(vmax == None):
        vmax = np.nanmax(pdata)

    if(ptype == 'trend'):
        colormap = plt.cm.bwr
    else:
        colormap = plt.cm.jet

    # Make copy of OMI_data array
    local_data  = np.copy(pdata)
    # Grid the data, fill in white space
    cyclic_data,cyclic_lons = add_cyclic_point(local_data,plon[0,:])
    plat,plon = np.meshgrid(plat[:,0],cyclic_lons)   
  
    # Mask any missing values
    mask_AI = np.ma.masked_where(cyclic_data < -998.9, cyclic_data)
    mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)

    #pax.gridlines()
    pax.coastlines(resolution='50m')
    mesh = pax.pcolormesh(plon, plat,\
            mask_AI.T,transform = datacrs,\
            cmap = colormap,vmin=vmin,vmax=vmax)
    pax.set_extent([-180,180,minlat,90],datacrs)
    pax.set_boundary(circle, transform=pax.transAxes)
    #cbar = plt.colorbar(mesh,ticks = np.arange(-200.,400.,5.),\
    cbar = plt.colorbar(mesh,\
        ax = pax, orientation='vertical',shrink = 0.8, extend = 'both')
    cbar.set_label(plabel,fontsize=colorbar_label_size,weight='bold')
    print("USING ",ptitle)
    pax.set_title(ptitle)


# Plot a single daily average
def plotCERES_Daily(CERES_data,day_idx,minlat = 60,save=False):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Make figure title
    title = 'CERES ' + CERES_data['param'] + ' ' + CERES_data['dates'][day_idx]

    # Mask the data
    local_data = np.copy(CERES_data['data'][day_idx,:,:])
    local_mask = np.ma.masked_where(local_data[:,:] < 0, local_data[:,:])

    # Make figure
    plt.close()
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(CERES_data['lon'], CERES_data['lat'],\
            local_mask,transform = datacrs,\
            cmap = colormap)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(0.0,400.,25.),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label=CERES_data['parm_name'])
    ax.set_title(title)

    if(save == True):   
        outname = 'ceres_grid_'+CERES_data['satellite'].lower()+'_'+\
                CERES_data['param']+'_'+CERES_data['dates'][day_idx]+'.png'
        plt.savefig(outname)
        print("Saved image",outname)
    else:
        plt.show()

# Plot a single monthly average
def plotCERES_Month(CERES_data, month_idx, pax = None, minlat = 65, \
        save = False):

    local_data = np.copy(CERES_data['data'][month_idx,:,:])
    mask_data = np.ma.masked_where((local_data < 0) | \
        (CERES_data['lat'] < minlat),local_data)

    # Make figure title
    title = 'CERES ' + CERES_data['param'] + ' ' + CERES_data['dates'][month_idx]

    if(pax is None):
        plt.close('all')
        fig1 = plt.figure()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)
 
        plotCERES_spatial(ax, CERES_data['lat'], CERES_data['lon'], \
            mask_data, 'climo', ptitle = title, plabel = 'W/m2', \
            vmin = None, vmax = None, colorbar_label_size = 14, \
            minlat = minlat)

        fig1.tight_layout()

        if(save == True):
            outname = 'ceres_grid_'+CERES_data['satellite'].lower()+'_'+\
                CERES_data['param']+'_'+CERES_data['dates'][month_idx]+'.png'
            plt.savefig(outname,dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

    else:
        plotCERES_spatial(pax, CERES_data['lat'], CERES_data['lon'], \
            mask_data, 'climo', ptitle = title, plabel = 'W/m2', \
            vmin = None, vmax = None, colorbar_label_size = 14, \
            minlat = minlat)

# Plot a monthly climatology 
def plotCERES_MonthClimo(CERES_data,month_idx,minlat = 60):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    month_obj = datetime(year = 1,month = int(CERES_data['dates'][month_idx][4:]),day = 1)
    str_month = month_obj.strftime('%B')

    # Make figure title
    title = 'CERES ' + CERES_data['param'] + ' ' + str_month + \
        ' Climatology\n'+str_month+'. ' + CERES_data['dates'][0][:4] \
        + ' - '+str_month+'. ' + CERES_data['dates'][-1][:4]

    # Make figure
    plt.close()
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(CERES_data['lon'], CERES_data['lat'],\
            CERES_data['MONTH_CLIMO'][month_idx,:,:],transform = datacrs,\
            cmap = colormap)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_boundary(circle, transform=ax.transAxes)
    #ax.set_xlim(-3430748.535086173,3430748.438879491)
    #ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(0.0,400.,25.),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.845,label=CERES_data['parm_name'])
    ax.set_title(title)

    plt.show()

# Designed to work with the netCDF data
def plotCERES_MonthTrend(CERES_data,month_idx=None,save=False,\
        trend_type='standard',season='',minlat=65.,return_trend=False, \
        pax = None):

    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='_thielSen'

    if(month_idx == None):
        month_adder = ''
        month_idx = 0
        index_jumper = 1
        do_month = False
        v_max = 30.
        v_min = -30.
    else:
        month_adder = '_month'
        if(CERES_data['season'] == 'sunlight'):
            index_jumper = 6
        else:   
            index_jumper = 12
        do_month = True
        v_max = 30.
        v_min = -30.

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(0,360,1.0)

    # --------------------------------------------------------------
    #
    # Use calcCERES_grid_trend to calculate the trends in the AI data
    #
    # --------------------------------------------------------------
    ceres_trends = calcCERES_grid_trend(CERES_data, month_idx, trend_type, \
        minlat)

    # --------------------------------------------------------------
    #
    # Plot the calculated trends on a figure
    #
    # --------------------------------------------------------------

    # Set up mapping variables 
    colormap = plt.cm.bwr

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(CERES_data['dates'][month_idx::index_jumper][0],'%Y%m')
    end_date   = datetime.strptime(CERES_data['dates'][month_idx::index_jumper][-1],'%Y%m')
    
    # Make figure title
    #date_month = datetime(year = 1,month = month_idx+1, day = 1).strftime('%B')
    month_string = ''
    if(do_month == True):
        month_string = start_date.strftime('%B') + ' '

    title = 'CERES ' + CERES_data['param'] + ' ' + month_string + 'Trends'\
        '\n'+start_date.strftime('%b. %Y') + ' - ' +\
        end_date.strftime('%b. %Y')

    # Call plotCERES_spatial to add the data to the figure

    if(pax is None):
        plt.close('all')
        fig1 = plt.figure()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

        plotCERES_spatial(ax, CERES_data['lat'], CERES_data['lon'], \
            ceres_trends, 'trend', ptitle = title, plabel = 'W/m2 per study period', \
            vmin = v_min, vmax = v_max, colorbar_label_size = 16, \
            minlat = minlat)

        if(save == True):
            month_adder = ''
            if(do_month == True):
                month_adder = '_' + start_date.strftime('%b') 
            out_name = 'ceres_' + CERES_data['param'] + '_trend'+ month_adder + '_' + \
                start_date.strftime('%Y%m') + '_' + end_date.strftime('%Y%m') + \
                '_min' + str(int(minlat)) + '.png'
            plt.savefig(out_name,dpi=300)
            print("Saved image",out_name)
        else:
            plt.show()
    else:
        plotCERES_spatial(pax, CERES_data['lat'], CERES_data['lon'], \
            ceres_trends, 'trend', ptitle = title, plabel = 'W/m2 per study period', \
            vmin = v_min, vmax = v_max, colorbar_label_size = 16, \
            minlat = minlat)

    if(return_trend == True):
        return ceres_trends

# Plot a daily monthly anomalies, which are daily data with that month's
# climatology subtracted 
def plotCERES_Daily_MonthAnomaly(CERES_data,CERES_daily_data,day_idx,minlat = 60,save=True):

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Determine which month to use based on the inputted day
    month_idx = int(CERES_daily_data['dates'][day_idx][4:6]) - 1

    # Calculate anomalies
    CERES_anom = CERES_daily_data['data'][day_idx,:,:] - CERES_data['MONTH_CLIMO'][month_idx,:,:]

    # Make figure title
    title = 'CERES ' + CERES_daily_data['param'] + ' Monthly Anomaly ' + CERES_daily_data['dates'][day_idx]

    # Make figure
    plt.close()
    fig1 = plt.figure(figsize = (8,8))
    ax = plt.axes(projection = mapcrs)
    ax.gridlines()
    ax.coastlines(resolution='50m')
    mesh = ax.pcolormesh(CERES_data['lon'], CERES_data['lat'],\
            CERES_anom,transform = datacrs,\
            cmap = colormap)
    ax.set_extent([-180,180,minlat,90],datacrs)
    ax.set_xlim(-3430748.535086173,3430748.438879491)
    ax.set_ylim(-3413488.8763307533,3443353.899053069)
    cbar = plt.colorbar(mesh,ticks = np.arange(-200,200.,25.),orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.905,label=CERES_data['parm_name'])
    ax.set_title(title)

    if(save == True):
        out_name = 'ceres_grid_'+CERES_daily_data['satellite'].lower() + '_' + CERES_daily_data['param'] + '_mnth_anom_'+\
            CERES_daily_data['dates'][day_idx]+'.png'       
        plt.savefig(out_name)
        print('Saved image '+out_name)
    else:
        plt.show()

# Makes a scatter plot comparing the trends in OMI and CERES
# ----------------------------------------------------------
def plot_compare_OMI_CERES_trends(OMI_data,CERES_data,month,minlat=65,\
        save=False):
    from OMILib import plotOMI_MonthTrend
    OMI_trend = plotOMI_MonthTrend(OMI_data,month_idx=month,save=True,\
                trend_type='standard',season='',minlat=minlat,\
                return_trend=True)
    CERES_trend = plotCERES_MonthTrend(CERES_data,month_idx=month,save=True,\
                trend_type='standard',season='',minlat=minlat,\
                return_trend=True)
    # Convert the index to a string using datetime
    if(month != None):
        dt_obj = datetime.strptime(OMI_data['DATES'][month],"%Y%m")
        title = 'OMI AI / CERES ' + CERES_data['param'] + '\n'+ \
            dt_obj.strftime("%b") + " Trend Comparison"
        outname = 'omi_ceres_trend_comp_'+dt_obj.strftime("%b")+'_'+\
            OMI_data['VERSION']+'vCERES.png'
    else:
        title = 'OMI AI / CERES ' + CERES_data['param'] + ' Trend Comparison'
        outname = 'omi_ceres_trend_comp_'+\
            OMI_data1['VERSION']+'vCERES_min' + str(minlat) + '.png'

    # Flip the CERES data to convert the longitudes from 0 - 360 to -180 - 180
    # ------------------------------------------------------------------------
    local_lon = np.copy(CERES_data['lon'])
    local_lat = np.copy(CERES_data['lat'])
    over_180  = np.where(CERES_data['lon'][0,:] >= 179.9999)
    under_180 = np.where(CERES_data['lon'][0,:] < 179.9999)

    for ii in range(local_lon.shape[0]):
        local_lon[ii,:] = np.concatenate([local_lon[ii,:][over_180] - 360.,\
            local_lon[ii,:][under_180]])
        local_lat[ii,:] = np.concatenate([local_lat[ii,:][over_180],\
            local_lat[ii,:][under_180]])
        CERES_trend[ii,:] = np.concatenate([CERES_trend[ii,:][over_180],\
            CERES_trend[ii,:][under_180]])
   
    print("Before removal, ",CERES_trend.shape)
 
    # First, mask any OMI data that are outside the bounds of the CERES data
    # Assume that the OMI data are larger than the CERES data
    # ----------------------------------------------------------------------
    OMI_trend[np.where(np.isnan(OMI_trend))] = -999.
    where_matching = np.where((OMI_data['LAT'][:,0] >= \
        (np.min(local_lat) - 0.5)) & (OMI_data['LAT'][:,0] <= \
        (np.max(local_lat) - 0.5)))[0]
    OMI_trend = OMI_trend[where_matching, :]

    mask_trend1 = np.array(OMI_trend[  (OMI_trend != 0) & (CERES_trend != 0) \
        & (OMI_trend != -999.) & (CERES_trend != -999.)])
    mask_trend2 = np.array(CERES_trend[(OMI_trend != 0) & (CERES_trend != 0) \
        & (OMI_trend != -999.) & (CERES_trend != -999.)])

    print(mask_trend1.shape, mask_trend2.shape)

    print("Pearson:  ",pearsonr(mask_trend1,mask_trend2))
    print("Spearman: ",spearmanr(mask_trend1,mask_trend2))

    xy = np.vstack([mask_trend1,mask_trend2])
    z = stats.gaussian_kde(xy)(xy)

    # Plot a somewhat-robust best fit line using Huber Regression
    # -----------------------------------------------------------
    x_scaler,y_scaler = StandardScaler(), StandardScaler()
    x_train = x_scaler.fit_transform(mask_trend1[...,None])
    y_train = y_scaler.fit_transform(mask_trend2[...,None])

    model = HuberRegressor(epsilon=1)
    model.fit(x_train,y_train)
    
    test_x = np.array([np.min(mask_trend1),np.max(mask_trend1)])
    predictions = y_scaler.inverse_transform(\
        model.predict(x_scaler.transform(test_x[...,None])))

    plt.close('all')
    fig1 = plt.figure()
    plt.scatter(mask_trend1,mask_trend2,c=z,s=8)
    plt.plot(test_x,predictions,color='tab:green',linestyle='--',label='Huber Fit')
    # Plot an unrobust fit line using linear regression
    # -------------------------------------------------
    plt.plot(np.unique(mask_trend1),np.poly1d(np.polyfit(mask_trend1,\
        mask_trend2,1))(np.unique(mask_trend1)),color='tab:orange',\
        linestyle='--',label='Polyfit Fit')

    ##if((month == 0) | (month == 1)):
    ##    plt.xlim(-0.5,0.3)
    ##    plt.ylim(-0.5,0.3)
    ##elif((month == 2)):
    ##    plt.xlim(-0.6,0.5)
    ##    plt.ylim(-0.6,0.5)
    ##elif((month == 3)):
    ##    plt.xlim(-0.5,0.5)
    ##    plt.ylim(-0.5,0.5)
    ##elif((month == 4)):
    ##    plt.xlim(-0.5,0.7)
    ##    plt.ylim(-0.5,0.7)
    ##elif((month == 5)):
    ##    plt.xlim(-0.5,0.5)
    ##    plt.ylim(-0.5,0.5)
    ##else:
    ##    plt.xlim(-0.3,0.3)
    ##    plt.ylim(-0.3,0.3)
    plt.legend()
    plt.xlabel(OMI_data['VERSION'])
    plt.ylabel('CERES')
    plt.title(title + '\nMinlat = ' + str(minlat))
    if(save == True):
        plt.savefig(outname)
    else:
        plt.show()
    #return mask_trend1,mask_trend2

def plot_compare_OMI_CERES_hrly(OMI_date,CERES_date,minlat=65,max_AI = -200.,\
        omi_dtype = 'control', only_sea_ice = False, skiprows = None, save=False):
#def plot_compare_OMI_CERES_hrly(OMI_hrly,CERES_hrly,minlat=65,save=False):

    if('/home/bsorenson/Research/OMI' not in sys.path):
        sys.path.append('/home/bsorenson/Research/OMI')
    from OMILib import readOMI_swath_shawn, readOMI_swath_hdf, \
        plotOMI_single_swath

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    plt.close('all')
    fig1 = plt.figure(figsize = (6,6))
    mapcrs = ccrs.NorthPolarStereo()
    ax0 = fig1.add_subplot(1,1,1, projection = mapcrs)

    # Set up mapping variables 
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.jet
    if(minlat < 45):
        mapcrs = ccrs.Miller()
    else:
        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)

    # Step 1: read in OMI and CERES data
    # ----------------------------------
    #OMI_hrly = readOMI_single_swath(OMI_date,60,only_sea_ice = only_sea_ice)
    if(omi_dtype == 'shawn'):
        OMI_hrly  = readOMI_swath_shawn(OMI_date, latmin = minlat,\
            skiprows = skiprows)
    else:
        if(omi_dtype == 'jz'):
            omi_dtype = 'JZ'
        OMI_hrly  = readOMI_swath_hdf(OMI_date, omi_dtype, \
            only_sea_ice = only_sea_ice, latmin = minlat, \
            skiprows = skiprows)

    CERES_hrly = readgridCERES_hrly(CERES_date,'SWF')

    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    plt.close('all')
    #fig, axs = plt.subplots(1,3,subplot_kw={'projection': mapcrs})
    fig = plt.figure(figsize=(17,5))
    ax0 = fig.add_subplot(1,3,2,projection = mapcrs)
    ax1 = fig.add_subplot(1,3,3,projection = mapcrs)
    ax2 = fig.add_subplot(1,3,1)

    # Step 3: Plot OMI and CERES data in first two panels
    # ---------------------------------------------------

    # -------------------------------------------------------
    # Use the single-swath plotting function to plot OMI data
    # -------------------------------------------------------
    plotOMI_single_swath(ax0, OMI_hrly, title = omi_dtype.title(), \
        circle_bound = True, gridlines = False)

    ax0.set_extent([-180,180,minlat,90], datacrs)


    local_data = np.copy(CERES_hrly['data'])

    # -------------------------------------------------------
    # Use the single-swath plotting function to plot OMI data
    # -------------------------------------------------------
    plot_lat, plot_lon = np.meshgrid(CERES_hrly['lat'],CERES_hrly['lon'])
    mask_flux = np.ma.masked_where(((CERES_hrly['counts'] == 0) | \
        (CERES_hrly['lat'] < minlat)), local_data)

    ax1.gridlines(ylocs = np.arange(minlat,90,5))
    ax1.coastlines(resolution = '50m')
    ax1.set_title('CERES ' + CERES_hrly['param'] + ' ' + CERES_hrly['date'])
    #plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
    mesh1 = ax1.pcolormesh(plot_lon, plot_lat,mask_flux,transform = datacrs,cmap = colormap,\
            vmin = min_dict[CERES_hrly['param']], \
            vmax = max_dict[CERES_hrly['param']])
    ax1.set_extent([-180,180,minlat,90],datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
            #vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])
    cbar = plt.colorbar(mesh1,ax=ax1,orientation='vertical',pad=0.02,\
        aspect=50,shrink = 1.000,label=CERES_hrly['param'] + ' [W/m2]')
    
    # Step 4: Plot scatter OMI and CERES comparison in third panel
    # ------------------------------------------------------------

    # Make gridded lat and lon arrays to use for coloring the plot
    glat, glon = np.meshgrid(OMI_hrly['LAT'],OMI_hrly['LON'])
    
    # Mask any empty boxes
    mask_AI = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
        (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),OMI_hrly['AI'])
    mask_flux = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
        (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),CERES_hrly['data'])
    mask_sza  = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
        (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),OMI_hrly['SZA'])
    #mask_lat  = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
    #    (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),glat)
    ##!## Convert the index to a string using datetime
    ##!#if(month != None):
    ##!#    dt_obj = datetime.strptime(OMI_data['DATES'][month],"%Y%m")
    ##!#    title = 'OMI AI / CERES ' + CERES_data['param'] + '\n'+ \
    ##!#        dt_obj.strftime("%b") + " Trend Comparison"
    ##!#    outname = 'omi_ceres_trend_comp_'+dt_obj.strftime("%b")+'_'+\
    ##!#        OMI_data['VERSION']+'vCERES.png'
    ##!#else:
    ##!#    title = 'OMI AI / CERES ' + CERES_data['param'] + ' Trend Comparison'
    ##!#    outname = 'omi_ceres_trend_comp_'+\
    ##!#        OMI_data1['VERSION']+'vCERES.png'


    #print("Pearson:  ",pearsonr(mask_AI,mask_flux))
    #print("Spearman: ",spearmanr(mask_AI,mask_flux))

    ##!#xy = np.vstack([mask_AI,mask_flux])
    ##!#z = stats.gaussian_kde(xy)(xy)

    ##!## Plot a somewhat-robust best fit line using Huber Regression
    ##!## -----------------------------------------------------------
    ##!#x_scaler,y_scaler = StandardScaler(), StandardScaler()
    ##!#x_train = x_scaler.fit_transform(mask_AI[...,None])
    ##!#y_train = y_scaler.fit_transform(mask_flux[...,None])

    ##!#model = HuberRegressor(epsilon=1)
    ##!#model.fit(x_train,y_train)
    ##!#
    ##!#test_x = np.array([np.min(mask_AI),np.max(mask_AI)])
    ##!#predictions = y_scaler.inverse_transform(\
    ##!#    model.predict(x_scaler.transform(test_x[...,None])))

    ##!## One to one line stuff
    ##!#xs = np.arange(np.min(mask_AI),np.max(mask_AI),0.1)

    #plt.close('all')
    #fig1 = plt.figure()
    #scat = ax2.scatter(mask_AI,mask_flux,c=mask_lat,s=8)
    scat = ax2.scatter(mask_AI,mask_flux,c=mask_sza,s=8)
    cbar = plt.colorbar(scat,ax=ax2,orientation='vertical',\
        label='Solar Zenith Angle',pad=0.02,aspect=50)
    ax2.set_title(OMI_hrly['date'])
    #plt.title(OMI_hrly['date'])
    #plt.scatter(mask_AI,mask_flux,c=z,s=8)
    #plt.plot(test_x,predictions,color='tab:green',linestyle='--',label='Huber Fit')
    #plt.plot(test_x,predictions,color='tab:green',linestyle='--',label='Huber Fit')
    # Plot an unrobust fit line using linear regression
    # -------------------------------------------------
    #plt.plot(np.unique(mask_AI),np.poly1d(np.polyfit(mask_AI,\
    #    mask_flux,1))(np.unique(mask_AI)),color='tab:orange',\
    #    linestyle='--',label='Polyfit Fit')
    # Plot a one-to-one line
    #plt.plot(xs,xs,label='1-1',color='tab:red')
    ##if((month == 0) | (month == 1)):
    ##    plt.xlim(-0.5,0.3)
    ##    plt.ylim(-0.5,0.3)
    ##elif((month == 2)):
    ##    plt.xlim(-0.6,0.5)
    ##    plt.ylim(-0.6,0.5)
    ##elif((month == 3)):
    ##    plt.xlim(-0.5,0.5)
    ##    plt.ylim(-0.5,0.5)
    ##elif((month == 4)):
    ##    plt.xlim(-0.5,0.7)
    ##    plt.ylim(-0.5,0.7)
    ##elif((month == 5)):
    ##    plt.xlim(-0.5,0.5)
    ##    plt.ylim(-0.5,0.5)
    ##else:
    ##    plt.xlim(-0.3,0.3)
    ##    plt.ylim(-0.3,0.3)
    #axs[2].legend()
    ax2.set_xlabel('OMI AI')
    ax2.set_ylabel('CERES ' + CERES_hrly['param'])
    #plt.subplots_adjust(wspace=0.1, hspace=0.10)
    plt.tight_layout()
    #plt.title(title)
    outname = 'omi_ceres_compare_'+OMI_date+'.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()
    return mask_AI,mask_flux

#def plot_compare_OMI_CERES_hrly_grid(OMI_date,CERES_date,minlat=65,max_AI = -200.,\
def plot_compare_OMI_CERES_hrly_grid(date_str,minlat=65,max_AI = -200.,\
        omi_dtype = 'control', only_sea_ice = False, only_ice = False, \
        no_ice = False,  skiprows = None, save=False, composite = False, \
        zoom = False, lat_circles = None):
#def plot_compare_OMI_CERES_hrly(OMI_hrly,CERES_hrly,minlat=65,save=False):

    if('/home/bsorenson/Research/OMI' not in sys.path):
        sys.path.append('/home/bsorenson/Research/OMI')
    if('/home/bsorenson/Research/MODIS' not in sys.path):
        sys.path.append('/home/bsorenson/Research/MODIS/obs_smoke_forcing')
    from OMILib import readOMI_swath_shawn, readOMI_swath_hdf, \
        plotOMI_single_swath
    from MODISLib import read_true_color

    # ------------------------------------------------------------------------
    #
    # Extract the OMI and CERES dates from the event date
    #
    # ------------------------------------------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    OMI_date = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                   dt_date_str.strftime('%H%M')]['omi']
    OMI_date = OMI_date.split('_')[-2][:4] + OMI_date.split('_')[-2][5:9] + \
        OMI_date.split('_')[-2][10:14]
    CERES_date = dt_date_str.strftime('%Y%m%d') + \
        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
        dt_date_str.strftime('%H%M')]['ceres_time']

    # ------------------------------------------------------------------------
    #
    # Set up mapping variables 
    #
    # ------------------------------------------------------------------------
    datacrs = ccrs.PlateCarree() 
    colormap = plt.cm.viridis
    circle_bound = True
    zoom_add = ''
    if(zoom):
        zoom_add = '_zoom'
        circle_bound = False
        mapcrs = ccrs.NorthPolarStereo(central_longitude = -40)
    else:
        mapcrs = ccrs.NorthPolarStereo()

    # ------------------------------------------------------------------------
    #
    # Step 1: read in OMI, CERES, and MODIS data
    #
    # ------------------------------------------------------------------------
    #OMI_hrly = readOMI_single_swath(OMI_date,60,only_sea_ice = only_sea_ice)
    if(omi_dtype == 'shawn'):
        print("Reading Shawn data")
        OMI_base  = readOMI_swath_shawn(OMI_date, latmin = minlat)
    else:
        OMI_base  = readOMI_swath_hdf(OMI_date, omi_dtype, \
            only_sea_ice = only_sea_ice, only_ice = only_ice, \
            no_ice = no_ice, latmin = minlat, \
            skiprows = skiprows)

    CERES_hrly = readgridCERES_hrly_grid(CERES_date, 'SWF', minlat = 55.)

    ##!## Read the true color data
    ##!## ------------------------
    ##!#var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
    ##!#    composite=composite)


    # ------------------------------------------------------------------------
    #
    # Step 2: Set up figure to have 3 panels
    #
    # ------------------------------------------------------------------------
    plt.close('all')
    #fig, axs = plt.subplots(1,3,subplot_kw={'projection': mapcrs})
    fig = plt.figure(figsize=(10,8))
    mapcrs = init_proj(date_str)
    ax0 = fig.add_subplot(2,2,1,projection = mapcrs)
    ax1 = fig.add_subplot(2,2,2,projection = mapcrs)
    #ax3 = fig.add_subplot(2,2,3,projection = crs1)
    ax2 = fig.add_subplot(2,2,4)

    ##!## Plot the true-color data for the previous date
    ##!## ----------------------------------------------
    ##!#ax3.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
    ##!#    var1.y[-1], var1.y[0]), origin='upper')

    ##!## Zoom in the figure if desired
    ##!## -----------------------------
    ##!#if(zoom):
    ##!#    ax3.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
    ##!#        lat_lims1[1]],crs = datacrs)

    # ------------------------------------------------------------------------
    #
    # Step 3: Plot OMI and CERES data in first two panels
    #
    # ------------------------------------------------------------------------
    lat_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                dt_date_str.strftime('%H%M')]['Lat']
    lon_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                dt_date_str.strftime('%H%M')]['Lon']
    if(zoom):
        OMI_base['UVAI'] = np.ma.masked_where(~((OMI_base['LAT'] >= lat_lims[0]) & \
                                                (OMI_base['LAT'] <  lat_lims[1]) & \
                                                (OMI_base['LON'] >= lon_lims[0]) & \
                                                (OMI_base['LON'] <  lon_lims[1])), OMI_base['UVAI'])
        #OMI_base['UVAI'] = np.ma.masked_where(OMI_base['UVAI'] > 1, OMI_base['UVAI'])
        CERES_hrly['data'] = np.ma.masked_where(~((CERES_hrly['lat'] >= lat_lims[0]) & \
                                                  (CERES_hrly['lat'] <  lat_lims[1]) & \
                                                  (CERES_hrly['lon'] >= lon_lims[0]) & \
                                                  (CERES_hrly['lon'] <  lon_lims[1])), CERES_hrly['data'])

    # Extract OMI and CERES values at a desired point
    # -----------------------------------------------
    tlat = 85.1155
    tlon = -68.2791
    c_idx = nearest_gridpoint(tlat, tlon,\
        CERES_hrly['lat'], CERES_hrly['lon'])
    if(len(c_idx[0]) > 1):
        c_idx = (np.array([c_idx[0][0]])), (np.array([c_idx[1][0]]))
    tflux = CERES_hrly['data'][c_idx]
    tfluxlat = CERES_hrly['lat'][c_idx]
    tfluxlon = CERES_hrly['lon'][c_idx]

    o_idx = nearest_gridpoint(tlat, tlon,\
        OMI_base['LAT'], OMI_base['LON'])
    if(len(o_idx[0]) > 1):
        print('o_idx = ',o_idx)
        o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
    tai    = OMI_base['UVAI'][o_idx]
    tailat = OMI_base['LAT'][o_idx]
    tailon = OMI_base['LON'][o_idx]

    ax0.plot(tailon, tailat,
             color='k', linewidth=2, marker='o',
             transform=datacrs)
    ax0.plot(tfluxlon, tfluxlat,
             color='r', linewidth=2, marker='o',
             transform=datacrs)
    ax0.plot(tlon, tlat,
             color='w', linewidth=2, marker='o',
             transform=datacrs)
    ax1.plot(tailon, tailat,
             color='k', linewidth=2, marker='o',
             transform=datacrs)
    ax1.plot(tfluxlon, tfluxlat,
             color='r', linewidth=2, marker='o',
             transform=datacrs)
    ax1.plot(tlon, tlat,
             color='w', linewidth=2, marker='o',
             transform=datacrs)

    print("Selected OMI and CERES at ", tlat, ' ', tlon,':')
    print(' OMI   - ',tai)
    print('     lat - ', tailat)
    print('     lon - ', tailon)
    print(' CERES - ',tflux)
    print('     lat - ', tfluxlat)
    print('     lon - ', tfluxlon)

    # Find averages of the data within the ergion


    # Use the single-swath plotting function to plot OMI data
    # -------------------------------------------------------
    plotOMI_single_swath(ax0, OMI_base, title = 'OMI ' + omi_dtype + \
        ' ' + OMI_date, circle_bound = circle_bound, gridlines = False)

    # Use the single-swath plotting function to plot CERES data
    # ---------------------------------------------------------
    plotCERES_hrly(ax1, CERES_hrly, minlat = minlat, \
        vmin = None, vmax = None, title = '', label = '', \
        circle_bound = circle_bound, gridlines = False, grid_data = True)

    if(zoom):
        lat_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                    dt_date_str.strftime('%H%M')]['Lat']
        lon_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                    dt_date_str.strftime('%H%M')]['Lon']

        ax0.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], datacrs)
        ax1.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], datacrs)
    else: 
        ax0.set_extent([-180,180,minlat,90], datacrs)
        ax1.set_extent([-180,180,minlat,90], datacrs)

    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax0, lat_circles) 
        plot_lat_circles(ax1, lat_circles) 

    # ------------------------------------------------------------------------
    #
    # Step 4: Colocate the OMI and CERES data
    #
    # ------------------------------------------------------------------------

    # First, get rid of any data outside the study region
    # ---------------------------------------------------
    print('Before zoom - ', OMI_base['UVAI'].shape)
    if(zoom):
        zoom_omi_LAT  = OMI_base['LAT'][( (OMI_base['LAT'] >= lat_lims[0]) & \
                                          (OMI_base['LAT'] <  lat_lims[1] ) & \
                                          (OMI_base['LON'] >= lon_lims[0]) & \
                                          (OMI_base['LON'] <  lon_lims[1] ))]
        zoom_omi_LON  = OMI_base['LON'][( (OMI_base['LAT'] >= lat_lims[0]) & \
                                          (OMI_base['LAT'] <  lat_lims[1] ) & \
                                          (OMI_base['LON'] >= lon_lims[0]) & \
                                          (OMI_base['LON'] <  lon_lims[1] ))]
        zoom_omi_UVAI = OMI_base['UVAI'][((OMI_base['LAT'] >= lat_lims[0]) & \
                                          (OMI_base['LAT'] <  lat_lims[1] ) & \
                                          (OMI_base['LON'] >= lon_lims[0]) & \
                                          (OMI_base['LON'] <  lon_lims[1] ))]
        ##!#zoom_ceres_LAT  = CERES_hrly['lat'][( (CERES_hrly['lat'] >= lat_lims[0]) & \
        ##!#                                      (CERES_hrly['lat'] <  lat_lims[1] ) & \
        ##!#                                      (CERES_hrly['lon'] >= lon_lims[0]) & \
        ##!#                                      (CERES_hrly['lon'] <  lon_lims[1] ))]
        ##!#zoom_ceres_LON  = CERES_hrly['lon'][( (CERES_hrly['lat'] >= lat_lims[0]) & \
        ##!#                                      (CERES_hrly['lat'] <  lat_lims[1] ) & \
        ##!#                                      (CERES_hrly['lon'] >= lon_lims[0]) & \
        ##!#                                      (CERES_hrly['lon'] <  lon_lims[1] ))]
        ##!#zoom_ceres_UVAI = CERES_hrly['data'][((CERES_hrly['lat'] >= lat_lims[0]) & \
        ##!#                                      (CERES_hrly['lat'] <  lat_lims[1] ) & \
        ##!#                                      (CERES_hrly['lon'] >= lon_lims[0]) & \
        ##!#                                      (CERES_hrly['lon'] <  lon_lims[1] ))]
        ##!#zoom_ceres_UVAI = CERES_hrly['data'][((CERES_hrly['lat'] >= lat_lims[0]) & \
        ##!#                                      (CERES_hrly['lat'] <  lat_lims[1] ) & \
        ##!#                                      (CERES_hrly['lon'] >= lon_lims[0]) & \
        ##!#                                      (CERES_hrly['lon'] <  lon_lims[1] ))]
    else:
        zoom_omi_LAT  = OMI_base['LAT']
        zoom_omi_LON  = OMI_base['LON']
        zoom_omi_UVAI = OMI_base['UVAI']
        ##!#zoom_ceres_LAT  = CERES_hrly['lat']
        ##!#zoom_ceres_LON  = CERES_hrly['lon']
        ##!#zoom_ceres_UVAI = CERES_hrly['data']
                                             

    # First, get rid of any masked values in the OMI data
    # ---------------------------------------------------
    mask_LAT  = zoom_omi_LAT[ ~zoom_omi_UVAI.mask]
    mask_LON  = zoom_omi_LON[ ~zoom_omi_UVAI.mask]
    mask_UVAI = zoom_omi_UVAI[~zoom_omi_UVAI.mask]
    mask_UVAI = np.ma.masked_where(mask_UVAI < -2, mask_UVAI)
    ##!#mask_LAT  = OMI_base['LAT'][~OMI_base['UVAI'].mask]

    ##!#mask_LON  = OMI_base['LON'][~OMI_base['UVAI'].mask]
    ##!#mask_UVAI = OMI_base['UVAI'][~OMI_base['UVAI'].mask]
    print('After zoom - ', mask_UVAI.shape)

    ceres_match_lat  = np.full(mask_UVAI.shape,np.nan)
    ceres_match_lon  = np.full(mask_UVAI.shape,np.nan)
    ceres_match_flux = np.full(mask_UVAI.shape,np.nan)
    ceres_match_vza  = np.full(mask_UVAI.shape,np.nan)

    # Loop over the OMI data and find the closest CERES
    # gridpoint
    # -------------------------------------------------
    for ii in range(ceres_match_flux.shape[0]):
        # Find the gridpoint in the gridded lat/lon data that 
        # corresponds to the station at slat and slon
        # ---------------------------------------------------- 
        o_idx = nearest_gridpoint(mask_LAT[ii], mask_LON[ii],\
            CERES_hrly['lat'], CERES_hrly['lon'])

        if(len(o_idx[0]) > 1):
            o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
        ceres_match_flux[ii] = CERES_hrly['data'][o_idx]
        ceres_match_lat[ii] = CERES_hrly['lat'][o_idx] 
        ceres_match_lon[ii] = CERES_hrly['lon'][o_idx] 
        ceres_match_vza[ii] = CERES_hrly['vza'][o_idx] 
        dist = ((ceres_match_lat[ii] - mask_LAT[ii])**2. + \
            (ceres_match_lon[ii] - mask_LON[ii])**2.)**0.5
        if( dist > 1.00):
        #    print(mask_LAT[ii], mask_LON[ii], mask_UVAI[ii], CERES_hrly['data'][o_idx], dist, ' mismatch')
            ceres_match_flux[ii] = np.nan
            ceres_match_lat[ii]  = np.nan
            ceres_match_lon[ii]  = np.nan
            ceres_match_vza[ii]  = np.nan
        #else:    
        #    print(mask_LAT[ii], mask_LON[ii], mask_UVAI[ii], CERES_hrly['data'][o_idx], dist)
            #print(mask_LAT[ii], mask_LON[ii], CERES_hrly['lat'][o_idx], CERES_hrly['lon'][o_idx])


    # ------------------------------------------------------------------------
    #
    # Step 5: Plot scatter OMI and CERES comparison in third panel
    #
    # ------------------------------------------------------------------------

    # Mask any empty boxes
    print('Before scatter', np.ma.masked_invalid(ceres_match_flux).compressed().shape)
    scat = ax2.scatter(mask_UVAI,ceres_match_flux,s=8)
    #scat = ax2.scatter(mask_AI,mask_flux,c=mask_sza,s=8)
    #cbar = plt.colorbar(scat,ax=ax2,orientation='vertical',\
    #    label='Viewing Zenith Angle',pad=0.02,aspect=50)
    ax2.set_title(OMI_base['date'])
    ax2.set_xlabel('OMI AI')
    ax2.set_ylabel('CERES ' + CERES_hrly['param'])
    fig.tight_layout()


    if(save == True):
        outname = 'omi_ceres_compare_'+OMI_date+'_TEST.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

    return OMI_base, CERES_hrly
