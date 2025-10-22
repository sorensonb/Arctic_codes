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
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
#import matplotlib.colors as color
import matplotlib.path as mpath
from matplotlib.colors import rgb2hex,Normalize
from matplotlib import cm
import matplotlib.dates as mdates
from matplotlib.cm import ScalarMappable
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os
from scipy.stats import pearsonr,spearmanr
from sklearn.linear_model import HuberRegressor
from sklearn.preprocessing import StandardScaler
from scipy.signal import find_peaks
import h5py
# The commands module was discontinued for Python 3, so if the user
# is using python 2, import commands instead

home_dir = os.environ['HOME']

sys.path.append(home_dir)
from python_lib import circle, plot_trend_line, nearest_gridpoint, \
    aerosol_event_dict, init_proj, plot_lat_circles, plot_figure_text, \
    plot_subplot_label, add_gridlines

if(home_dir + '/Research/NSIDC' not in sys.path):
    sys.path.append(home_dir + '/Research/NSIDC')
if(home_dir + '/Research/OMI' not in sys.path):
    sys.path.append(home_dir + '/Research/OMI')
if(home_dir + '/Research/NAAPS' not in sys.path):
    sys.path.append(home_dir + '/Research/NAAPS')
from NSIDCLib import *
from glob import glob
#from NAAPSLib import *
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
    'CLR': 100.,
    'CLD': 100.,
    'ALB': 1.
}

min_dict = {
    'SWF': 120.,
    'LWF': 300.,
    'TOTAL': 450.,
    'CLR': 0.,
    'CLD': 0.,
    'ALB': 0.
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
        base_path = home_dir + 'data/CERES/SSF_1Deg/monthly/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    else:
        base_path = home_dir + '/data/CERES/SSF_1Deg/monthly/Aqua/CERES_SSF1deg-Month_Aqua-MODIS_Ed4.1_Subset_'
    #base_path = '/data/CERES/SSF_1Deg/monthly/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    #base_path = '/home/bsorenson/data/CERES/SSF_1Deg/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    total_list = sorted(glob(base_path+'*.nc')) 

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

def readgridCERES_daily(date_str, end_str = None, satellite = 'Aqua', \
        minlat = 60.5):

    dt_begin_str = datetime.strptime(date_str, '%Y%m%d')
    dt_end_str = datetime.strptime(end_str, '%Y%m%d')
    dt_begin_noaa20 = datetime.strptime('20180401', '%Y%m%d')
    dt_begin_suominpp = datetime.strptime('20120401', '%Y%m%d')
    dt_end_suominpp = datetime.strptime('20190930', '%Y%m%d')

    if((satellite == 'Aqua') | (satellite == 'Terra') | \
        (satellite == 'SuomiNPP') | (satellite == 'NOAA20')):
        CERES_data = readgridCERES_daily_file(date_str, end_str = end_str, \
            satellite = satellite, minlat = minlat)
    else:
        
        found_NPP = False
        found_NOAA = False

        CERES_data1 = readgridCERES_daily_file(date_str, end_str = end_str, \
            satellite = 'Aqua', minlat = minlat)
        CERES_data2 = readgridCERES_daily_file(date_str, end_str = end_str, \
            satellite = 'Terra', minlat = minlat)
        if( (dt_begin_str > dt_begin_noaa20) ):
            found_NOAA = True
            CERES_data5 = readgridCERES_daily_file(date_str, end_str = end_str, \
                satellite = 'NOAA20', minlat = minlat)
        if( (dt_begin_str > dt_begin_suominpp) & (dt_end_str <= dt_end_suominpp)):
            found_NPP = True
            CERES_data3 = readgridCERES_daily_file(date_str, end_str = end_str, \
                satellite = 'SuomiNPP', minlat = minlat)
        CERES_data4 = readgridCERES_daily_file(date_str, end_str = end_str, \
            satellite = 'Aqua',minlat = minlat)

        pvars = ['alb_all', 'alb_clr', 'swf_all', 'swf_clr', 'lwf_all', \
            'lwf_clr', 'cld', 'ice_conc']

        over_180    = np.where(CERES_data1['lon'] < 0)
        under_180   = np.where(CERES_data1['lon'] >= 0)

        #for ii in range(local_lon.shape[0]):
        #    local_lon[ii,:] = np.concatenate([local_lon[ii,:][over_180] - 360.,\
        #        local_lon[ii,:][under_180]])

        for pvar in pvars:
            total_data = np.concatenate((CERES_data1[pvar], CERES_data2[pvar]))
            if(found_NPP):
                total_data = np.concatenate((total_data, CERES_data3[pvar]))
            if(found_NOAA):
                total_data = np.concatenate((total_data, CERES_data5[pvar]))
                #CERES_data3[pvar]))
                #CERES_data3[pvar], CERES_data5[pvar]))
            total_data = np.ma.masked_where(total_data < 0, total_data)
            #total_data = np.concatenate((total_data[over_180], \
            #    total_data[under_180]))

            CERES_data4[pvar] = total_data 

        #CERES_data4['lon'] = np.concatenate((CERES_data4['lon'][over_180], \
        #    CERES_data4['lon'][under_180]))
        CERES_data4['satellite'] = 'All' 

        return CERES_data4

# For now, assume that the user only wants to grab a single file of daily data
# User could cram multiple months into a single file, if wanted
#def readgridCERES_daily(start_date,end_date,param,satellite = 'Aqua',minlat=60.5,season='all'):
def readgridCERES_daily_file(date_str, end_str = None, satellite = 'Aqua', \
        minlat=60.5):
    CERES_data = {}
  
#    lat_ranges = np.arange(minlat,90.5,1.0)
#    lon_ranges = np.arange(0.5,360.5,1.0)

    # Grab all the files
    if(satellite == 'Terra'):
        base_path = home_dir + '/data/CERES/SSF_1Deg/daily/Terra/CERES_SSF1deg-Day_Terra-MODIS_Ed4.1_Subset_'
    elif(satellite == 'Aqua'):
        base_path = home_dir + '/data/CERES/SSF_1Deg/daily/Aqua/CERES_SSF1deg-Day_Aqua-MODIS_Ed4.1_Subset_'
    elif(satellite == 'SuomiNPP'):
        base_path = home_dir + '/data/CERES/SSF_1Deg/daily/SuomiNPP/CERES_SSF1deg-Day_NPP-VIIRS_Ed2A_Subset_'
    elif(satellite == 'NOAA20'):
        base_path = home_dir + '/data/CERES/SSF_1Deg/daily/NOAA20/CERES_SSF1deg-Day_NOAA20-VIIRS_Ed1B_Subset_'
    #base_path = '/data/CERES/SSF_1Deg/monthly/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    #base_path = '/home/bsorenson/data/CERES/SSF_1Deg/Terra/CERES_SSF1deg-Month_Terra-MODIS_Ed4A_Subset_'
    total_list = sorted(glob.glob(base_path+'*.nc'))

    # Loop over all files and find the ones that match with the desired times
    #start_date = datetime.strptime(str(start_date),'%Y%m%d')
    #end_date = datetime.strptime(str(end_date),'%Y%m%d')
    dt_date_str = datetime.strptime(date_str, '%Y%m%d')

    #final_list = []
    #for f in total_list:
    #    fdate = f.split('_')[-1][:6]
    #    fdate = datetime.strptime(str(fdate),'%Y%m')
    #    if((fdate>=start_date) & (fdate<=end_date)):
    #        final_list.append(f)
    # Loop over the list of current CERES data files and find the one that  
    # corresponds to the desired datetime
    # --------------------------------------------------------------------
    good_list = []
    for tfile in total_list:
        filename   = tfile.strip().split('/')[-1]
        begin_date = datetime.strptime(filename[-20:-12],"%Y%m%d")
        end_date   = datetime.strptime(filename[-11:-3], "%Y%m%d")
        if((dt_date_str >= begin_date) & (dt_date_str <= end_date)):
            print("Found matching file",tfile)
            good_list.append(tfile)
            #work_file = tfile

    # Open up the file
    data = Dataset(good_list[0],'r')

    # Convert the times to dt objects
    # -------------------------------
    base_date = datetime.strptime('20000301','%Y%m%d')
    dt_date_str = dt_date_str + timedelta(days = 0.5)
    times = np.array([base_date + timedelta(days = float(tdate)) for tdate in data['time']])
    if(end_str is None):
        time_idx = np.where(times == dt_date_str)[0]
    else:
        dt_end_str = datetime.strptime(end_str, '%Y%m%d') + timedelta(days = 0.5)
        time_idx = np.where( (times >= dt_date_str ) & (times <= dt_end_str))[0] 

    lat_idx = np.where(data['lat'][:] >= minlat)[0]

    newlon, newlat = np.meshgrid(data['lon'], data['lat'][lat_idx])
    dt_dates = times[time_idx]
    CERES_data['dt_dates'] = dt_dates
    CERES_data['alb_all']  = data['toa_alb_all_daily'][time_idx,lat_idx,:].squeeze()
    CERES_data['alb_clr']  = data['toa_alb_clr_daily'][time_idx,lat_idx,:].squeeze()
    CERES_data['swf_all']  = data['toa_sw_all_daily'][time_idx,lat_idx,:].squeeze()
    CERES_data['swf_clr']  = data['toa_sw_clr_daily'][time_idx,lat_idx,:].squeeze()
    CERES_data['lwf_all']  = data['toa_lw_all_daily'][time_idx,lat_idx,:].squeeze()
    CERES_data['lwf_clr']  = data['toa_lw_clr_daily'][time_idx,lat_idx,:].squeeze()
    CERES_data['cld']      = data['cldarea_total_day_daily'][time_idx,lat_idx,:].squeeze()
    CERES_data['ice_conc'] = data['aux_snow_daily'][time_idx,lat_idx,:].squeeze()
    CERES_data['lat'] = newlat
    CERES_data['lon'] = newlon
    #CERES_data['lat'] = data['lat'][lat_idx]
    #CERES_data['lon'] = data['lon'][:]
    #CERES_data['lon'][CERES_data['lon'] > 179.99] = -360. + \
    #    CERES_data['lon'][CERES_data['lon'] > 179.99]
    CERES_data['satellite'] = satellite

    #split_date = data.variables['time'].units.split()[2].split('-')
    #base_date = datetime(year = int(split_date[0]), month = int(split_date[1]),\
    #     day = int(split_date[2]))

    ## Loop over good files and insert data into dictionary
    #for ti in range(data.variables['time'].size):
    #    CERES_data['data'][ti,:,:] = data.variables[param][ti,lat_indices,:]
    #    new_date = base_date + relativedelta(days = data.variables['time'][ti])
    #    CERES_data['dates'].append(new_date.strftime("%Y%m%d"))
    #    ### print(data.variables[param][0,lat_indices,:])
    data.close() 
     
    #start_date = datetime.strptime(str(start_date),'%Y%m')
    #end_date = datetime.strptime(str(end_date),'%Y%m') +timedelta(days=31)
    #
    #tempdate = data.variables['time'].units
    #orig_date = datetime.strptime(tempdate.split()[2]+' '+tempdate.split()[3],'%Y-%m-%d %H:%M:%S')

    return CERES_data

# Data period is of format YYYYMMDDHH
# NOTE: SET modis_comp TO TRUE FOR THERMAL IR DIXIE STUDY
def readgridCERES_hrly_grid(data_dt,param,satellite = 'Aqua',minlat=60.0,\
        season='all', modis_comp = False):
    print(data_dt)
    lat_ranges = np.arange(minlat,90.0,0.25)
    lon_ranges = np.arange(-180.0,180.0,0.25)

    # Grab all the files
    if(satellite == 'Terra'):
        base_path = home_dir + '/data/CERES/SSF_Level2/Terra/'
    elif(satellite == 'Aqua'):
        base_path = home_dir + '/data/CERES/SSF_Level2/Aqua/'
        if(modis_comp):
            base_path = base_path + 'modis_comp/'
        #base_path = '/home/bsorenson/data/CERES/SSF_Level2/Aqua/modis_comp/'
    elif(satellite == 'SuomiNPP'):
        base_path = home_dir + '/data/CERES/SSF_Level2/SuomiNPP/'
    elif(satellite == 'NOAA20'):
        base_path = home_dir + '/data/CERES/SSF_Level2/NOAA20/'

    total_list = sorted(glob(base_path+'CERES_SSF_*.nc'))

    # Convert the desired dt to a datetime object to use for finding the file
    # -----------------------------------------------------------------------
    day = False
    if(len(data_dt) == 12):
        str_fmt = "%Y%m%d%H%M"
        hour_adder = 1
        hour_subtracter = 0
    elif(len(data_dt) == 10):
        str_fmt = "%Y%m%d%H"
        hour_adder = 1
        hour_subtracter = 0
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
    dt_data_end   = datetime.strptime(data_dt,str_fmt) \
        + relativedelta(hours = hour_adder)

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
    if(len(good_list) == 0):
        print("ERROR: could not find matching files. Returning")
        return -1

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
    if(satellite == 'NOAA20'):
        lat = data['Time_and_Position/instrument_fov_latitude'][:]
        lon = data['Time_and_Position/instrument_fov_longitude'][:]
        sza = data['Viewing_Angles/solar_zenith_angle'][:]
        vza = data['Viewing_Angles/view_zenith_angle'][:]
        azm = data['Viewing_Angles/relative_azimuth_angle'][:]
        alb = data['TOA_and_Surface_Fluxes/broadband_surface_albedo'][:]
        cld = 100. - data['Cloudy_Footprint_Area/layers_coverages'][:,0]
        lon[lon>179.99] = -360.+lon[lon>179.99]
        lwf = data['TOA_and_Surface_Fluxes/toa_longwave_flux'][:] 
        swf = data['TOA_and_Surface_Fluxes/toa_shortwave_flux'][:] 
    else:
        lat   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        lon   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        sza   = data.variables['CERES_solar_zenith_at_surface'][:]
        vza   = data.variables['CERES_viewing_zenith_at_surface'][:]
        azm   = data.variables['CERES_relative_azimuth_at_surface'][:]
        alb   = data.variables['CERES_broadband_surface_albedo'][:]
        cld   = 100. - data.variables['Clear_layer_overlap_percent_coverages'][:,0]
        #cld   = 100. - data.variables['Cloud_mask_clear_weak_percent_coverage'][:]
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
   
    data.close()
 
    # Loop over the values and rows
    print(swf.shape) 

    # Convert each pixel relative time to absolute
    total_times = np.array([base_date + relativedelta(days = ttime) \
        for ttime in time])

    ##!#plt.close('all') 
    ##!#fig1 = plt.figure()
    ##!#plt.plot(total_times, azm, label = 'azm')
    ##!#plt.plot(total_times, lat, label = 'lat')
    ##!#plt.plot(total_times, vza, label = 'vza')   
    ##!#plt.legend()
    ##!#plt.show()
   
    # Extract only the data in the time window 
    # ----------------------------------------
    test_time = total_times[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_ftim = time[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_swf  = swf[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lwf  = lwf[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lat  = lat[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lon  = lon[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_sza  = sza[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_vza  = vza[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_azm  = azm[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_alb  = alb[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_cld  = cld[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]

    # Determine where the LAT peaks are located and separated by 181
    # --------------------------------------------------------------
    lat_neg_peaks, _ = find_peaks(-test_vza, distance = 120)
    lat_peak_diffs = lat_neg_peaks[1:] - lat_neg_peaks[:-1]

    ##!#print('lat_peak_diffs = ', lat_peak_diffs)
    ##!#print('lat_neg_peak times = ', test_time[lat_neg_peaks])
    ##!#print('lat_neg_peak lats  = ', test_lat[lat_neg_peaks])

    ##!#print(lat_peak_diffs)
    ##!#plt.close('all') 
    ##!#fig1 = plt.figure()
    ##!#plt.plot(test_time, test_azm, label = 'azm')
    ##!#plt.plot(test_time, test_lat, label = 'lat')
    ##!#plt.plot(test_time, test_vza, label = 'vza')   
    ##!#plt.plot(test_time[lat_neg_peaks], test_azm[lat_neg_peaks], 'x')
    ##!#plt.plot(test_time[lat_neg_peaks], test_vza[lat_neg_peaks], 'x', color = 'k')
    ##!#plt.legend()
    ##!#plt.show()



    # Extract the data within the full cycles
    # ---------------------------------------    
    #keep_lat_peaks = lat_neg_peaks[np.where((lat_peak_diffs == 181) | \
    #    (lat_peak_diffs == 182))]
    keep_lat_peaks = lat_neg_peaks
    ##!#print('AA', lat_peak_diffs)
    ##!#print('keep_lat_peak times = ', test_time[keep_lat_peaks])
    ##!#print('keep_lat_peak lats  = ', test_lat[keep_lat_peaks])
    if(len(keep_lat_peaks) == 0):
        print("ERROR: no data within the desired time window")
        return

    begin_idx = keep_lat_peaks[0]
    end_idx   = keep_lat_peaks[-1]

    test_time = test_time[begin_idx: end_idx]
    test_ftim = test_ftim[begin_idx: end_idx]
    test_swf  = test_swf[begin_idx: end_idx]
    test_lwf  = test_lwf[begin_idx: end_idx]
    test_lat  = test_lat[begin_idx: end_idx]
    test_lon  = test_lon[begin_idx: end_idx]
    test_sza  = test_sza[begin_idx: end_idx]
    test_vza  = test_vza[begin_idx: end_idx]
    test_azm  = test_azm[begin_idx: end_idx]
    test_alb  = test_alb[begin_idx: end_idx]
    test_cld  = test_cld[begin_idx: end_idx]

    # Make grid arrays for the data
    # -----------------------------
    testmax = 181
    grid_swf  = np.full((300, testmax), np.nan)
    grid_lwf  = np.full((300, testmax), np.nan)
    grid_lat  = np.full((300, testmax), np.nan)
    grid_lon  = np.full((300, testmax), np.nan)
    grid_sza  = np.full((300, testmax), np.nan)
    grid_vza  = np.full((300, testmax), np.nan)
    grid_azm  = np.full((300, testmax), np.nan)
    grid_alb  = np.full((300, testmax), np.nan)
    grid_cld  = np.full((300, testmax), np.nan)
    grid_time = np.full((300, testmax), np.nan)

    # Loop over the data and insert into the grids
    for ii in range(len(keep_lat_peaks) - 1):
        idx_diff = keep_lat_peaks[ii+1] - keep_lat_peaks[ii]
        #print(idx_diff, test_time[keep_lat_peaks[ii]])
        if(idx_diff == testmax):
    ##!#        print(idx_diff, test_time[keep_lat_peaks[ii]], np.nanmean(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]), 'exact match')
            grid_swf[ii, :len(test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_lwf[ii, :len(test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            #tester = readgridCERES_hrly_grid(date_str[:10], param, satellite = 'Terra', minlat = 30)print(grid_lwf[ii, :len(test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])])
            grid_lat[ii, :len(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_lon[ii, :len(test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_sza[ii, :len(test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_vza[ii, :len(test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_azm[ii, :len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_alb[ii, :len(test_alb[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_alb[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_cld[ii, :len(test_cld[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_cld[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
            grid_time[ii,:len(test_ftim[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_ftim[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]
        ##!#elif(idx_diff > 181):
        ##!#    print(idx_diff, test_time[keep_lat_peaks[ii]], np.nanmean(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]), 'over estimate')
        ##!#    grid_swf[ii,:len(test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
        ##!#    grid_lwf[ii,:len(test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
        ##!#    grid_lat[ii,:len(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
        ##!#    grid_lon[ii,:len(test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
        ##!#    grid_sza[ii,:len(test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
        ##!#    grid_vza[ii,:len(test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
        ##!#    grid_azm[ii,:len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181])] = test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 181]
        elif(idx_diff < testmax):
        #else:
            #print(idx_diff, test_time[keep_lat_peaks[ii]], 'too small')
            # Find where the VZA in this row is closest to the VZA in the previous row
            # ------------------------------------------------------------------------

    ##!#        print(idx_diff, test_time[keep_lat_peaks[ii]], np.nanmean(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]), 'nothing')
            if(ii > 0):
                if(np.count_nonzero(~np.isnan(grid_lon[ii-1,:])) > 0):
                    if(keep_lat_peaks[ii] < len(test_vza)):
                        avg_vza = np.nanmean(grid_vza[:,:], axis=0)
                        fun_c = abs((test_vza[keep_lat_peaks[ii]]) - avg_vza[:30])
                        #fun_c = abs((test_vza[keep_lat_peaks[ii]]) - grid_vza[ii-1,:30])
                        m_idx = np.where(fun_c == np.min(fun_c))
                        ##!#print('last line vza = ',grid_vza[ii-1,:9])
                        ##!#print('avg line vza  = ',avg_vza[:9])
                        ##!##print('avg line vza  = ',np.nanmean(grid_vza[:,:], axis=0)[:9])
                        ##!#print('this line vza = ',test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + 9])
                        ##!#print('beginning lat this line = ', test_vza[keep_lat_peaks[ii]])
                        ##!#print('closest   lat last line = ', avg_vza[:30][m_idx])
                        ##!##print('closest   lat last line = ', grid_vza[ii-1,:9][m_idx])
                        ##!#print('index of closst last lat =', m_idx)
                        
                        #fun_c = np.maximum(np.abs(grid_lat - slat), \
                        #    np.abs(grid_lon - slon))
                        if(len(test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]) > 0):
                            grid_swf[ii,m_idx[0][0]:len(test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            grid_lwf[ii,m_idx[0][0]:len(test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            #print(grid_lwf[ii,m_idx[0][0]:len(test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])])
                            grid_lat[ii,m_idx[0][0]:len(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            grid_lon[ii,m_idx[0][0]:len(test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            grid_sza[ii,m_idx[0][0]:len(test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            grid_vza[ii,m_idx[0][0]:len(test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            grid_azm[ii,m_idx[0][0]:len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            grid_alb[ii,m_idx[0][0]:len(test_alb[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_alb[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            grid_cld[ii,m_idx[0][0]:len(test_cld[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_cld[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]
                            grid_time[ii,m_idx[0][0]:len(test_ftim[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]])] = test_ftim[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]][:-(int(m_idx[0][0]))]

        elif(idx_diff > testmax):
            #print(idx_diff, test_time[keep_lat_peaks[ii]], np.nanmean(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii+1]]), 'over estimate')
            grid_swf[ii,:len(test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_swf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_lwf[ii,:len(test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_lwf[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_lat[ii,:len(test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_lat[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_lon[ii,:len(test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_lon[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_sza[ii,:len(test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_sza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_vza[ii,:len(test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_vza[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_azm[ii,:len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_alb[ii,:len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_alb[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_cld[ii,:len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_cld[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
            grid_time[ii,:len(test_azm[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax])] = test_ftim[keep_lat_peaks[ii]:keep_lat_peaks[ii] + testmax]
        

    # Remove any rows in the grid arrays with any nans
    # ------------------------------------------------
    checkers = np.array([np.count_nonzero(~np.isnan(grid_lon[ii,:])) for ii in range(grid_lon.shape[0])])
    ##!#print(checkers) 

    ###grid_swf = grid_swf[np.where(checkers == 181)[0],:]
    ###grid_lwf = grid_lwf[np.where(checkers == 181)[0],:]
    ###grid_lat = grid_lat[np.where(checkers == 181)[0],:]
    ###grid_lon = grid_lon[np.where(checkers == 181)[0],:]
    ###grid_sza = grid_sza[np.where(checkers == 181)[0],:]
    ###grid_vza = grid_vza[np.where(checkers == 181)[0],:]
    ###grid_azm = grid_azm[np.where(checkers == 181)[0],:]
    grid_lat[np.isnan(grid_lat)] = -999.0
    grid_lon[np.isnan(grid_lon)] = -999.0

    grid_swf  = np.ma.masked_invalid(grid_swf)
    grid_lwf  = np.ma.masked_invalid(grid_lwf)
    grid_swf  = np.ma.masked_where((grid_swf > 5000), grid_swf)
    grid_lwf  = np.ma.masked_where((grid_lwf > 5000), grid_lwf)

    # Remove any rows that have missing geolocation data
    # --------------------------------------------------
    finders = np.array([-999.0 not in line for line in grid_lat])   

    grid_swf  = grid_swf[finders]
    grid_lwf  = grid_lwf[finders]
    grid_lat  = grid_lat[finders]
    grid_lon  = grid_lon[finders]
    grid_vza  = grid_vza[finders]
    grid_sza  = grid_sza[finders]
    grid_azm  = grid_azm[finders]
    grid_alb  = grid_alb[finders]
    grid_cld  = grid_cld[finders]
    grid_time = grid_time[finders]

    # Convert the day timedelta to actual datetimes
    # ---------------------------------------------
    grid_time_dt = np.array([[base_date + relativedelta(days = time) \
        for time in row] for row in grid_time])
 
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
    CERES_grid_hrly['alb']       = grid_alb
    CERES_grid_hrly['cld']       = grid_cld
    CERES_grid_hrly['time']      = grid_time
    CERES_grid_hrly['time_dt']   = grid_time_dt
    CERES_grid_hrly['base_date'] = '197001010000'
    CERES_grid_hrly['filename'] = good_list[0]

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

# Writes a MODIS channel dictionary to HDF5 for Fortran colocation
# NOTE: used for writing gridded SSF_Level2 to an HDF5 file for a
#       single hour.
# ----------------------------------------------------------------
def write_CERES_hrly_grid_to_HDF5(CERES_grid_hrly, save_path = './', \
        minlat = 65., remove_empty_scans = False):

    if(isinstance(CERES_grid_hrly, str)):
        CERES_grid_hrly =  readgridCERES_hrly_grid(CERES_grid_hrly, 'swf', \
            minlat=minlat)

    if(remove_empty_scans):
        mask_data = np.ma.masked_where(CERES_grid_hrly['lat'] < minlat, \
            CERES_grid_hrly['swf'])
        mask_dims = np.array([ (False in mask_data[ii,:].mask) for ii in \
            range(mask_data.shape[0])])
        keep_idxs = np.where(mask_dims == True)[0]
    else:
        keep_idxs = None 

    # Convert the filename object to datetime
    # ---------------------------------------
    file_date = CERES_grid_hrly['date']
    dt_date_str = datetime.strptime(file_date, '%Y%m%d%H')

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'ceres_subset_' + file_date + '.hdf5'
    dset = h5py.File(outfile,'w')
 
    dset.create_dataset('latitude',  data = CERES_grid_hrly['lat'][keep_idxs,:].squeeze())
    dset.create_dataset('longitude', data = CERES_grid_hrly['lon'][keep_idxs,:].squeeze())
    dset.create_dataset('swf',       data = CERES_grid_hrly['swf'][keep_idxs,:].squeeze())
    dset.create_dataset('lwf',       data = CERES_grid_hrly['lwf'][keep_idxs,:].squeeze())
    dset.create_dataset('cld',       data = CERES_grid_hrly['cld'][keep_idxs,:].squeeze())
    dset.create_dataset('alb',       data = CERES_grid_hrly['alb'][keep_idxs,:].squeeze())
    dset.create_dataset('time',      data = CERES_grid_hrly['time'][keep_idxs,:].squeeze())

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  

# For a given large CERES SSF_Level2 file containing, split out the data and 
# write each time period for each case. 
# NOTE: Used in the NAAPS/CERES comparisons
# ----------------------------------------------------------------------------
def write_CERES_L2_to_HDF5(case_date, satellite, save_path = './'):

    if(home_dir + '/Research/NAAPS/' not in sys.path):
        sys.path.append(home_dir + '/Research/NAAPS/')
    from NAAPSLib import event_dict

    # Use the case date to pick the time periods
    # ------------------------------------------
    before_start = event_dict[case_date]['before_start']
    before_end   = event_dict[case_date]['before_end']
    event_start  = event_dict[case_date]['start']
    event_end    = event_dict[case_date]['end']
    end_start    = event_dict[case_date]['end_start']
    end_end      = event_dict[case_date]['end_end']

    # Select and read the matching CERES file
    # ---------------------------------------

    data_dir = home_dir + '/data/CERES/SSF_Level2/' + satellite + '/'

    ## Grab all the files
    #if(satellite == 'Terra'):
    #    base_path = home_dir + '/data/CERES/SSF_Level2/Terra/'
    #elif(satellite == 'Aqua'):
    #    base_path = home_dir + '/data/CERES/SSF_Level2/Aqua/'
    #elif(satellite == 'SuomiNPP'):
    #    base_path = home_dir + '/data/CERES/SSF_Level2/SuomiNPP/'
    #elif(satellite == 'NOAA20'):
    #    base_path = home_dir + '/data/CERES/SSF_Level2/NOAA20/'

    total_list = sorted(glob.glob(data_dir + 'CERES_SSF_*.nc'))

    # Set up the before and after datetime objects
    # --------------------------------------------
    dt_begin_start = datetime.strptime(before_start, '%Y%m%d')
    dt_begin_end   = datetime.strptime(before_end,   '%Y%m%d')
    dt_event_start = datetime.strptime(event_start, '%Y%m%d%H')
    dt_event_end   = datetime.strptime(event_end,   '%Y%m%d%H')
    dt_end_start   = datetime.strptime(end_start, '%Y%m%d')
    dt_end_end     = datetime.strptime(end_end,   '%Y%m%d')

    # Loop over the list of current CERES data files and find the one that  
    # corresponds to the desired datetime
    # --------------------------------------------------------------------
    good_list = []
    for tfile in total_list:
        split_file = tfile.strip().split('/')[-1].split('_')
        begin_date = datetime.strptime(split_file[-1][:10],"%Y%m%d%H")
        end_date   = datetime.strptime(split_file[-1][11:21],"%Y%m%d%H")
        if((dt_begin_start >= begin_date) & (dt_end_end <= end_date)): # | \
           #(dt_data_begin == end_date) | \
           #((day == True) & (begin_date >= dt_data_begin) & \
           #(end_date <= dt_data_end))):
            print("Found matching file",tfile)
            good_list.append(tfile)
            #work_file = tfile
    if(len(good_list) == 0):
        print("ERROR: could not find matching files. Returning")
        return -1

    # Read in the CERES data from the file
    # ------------------------------------
    base_date = datetime(year=1970,month=1,day=1)

    # Assume that only one file will match up
    # --------------------------------------- 
    #for fileI in range(len(good_list)):
    data = Dataset(good_list[0],'r')
    time      = data.variables['time'][:]
    lat       = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
    lon       = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
    sza       = data.variables['CERES_solar_zenith_at_surface'][:]
    vza       = data.variables['CERES_viewing_zenith_at_surface'][:]
    azm       = data.variables['CERES_relative_azimuth_at_surface'][:]
    alb       = data.variables['CERES_broadband_surface_albedo'][:]
    lwf       = data.variables['CERES_LW_TOA_flux___upwards'][:]
    swf       = data.variables['CERES_SW_TOA_flux___upwards'][:]
    clr_str   = data.variables['Cloud_mask_clear_strong_percent_coverage'][:]
    clr_wek   = data.variables['Cloud_mask_clear_weak_percent_coverage'][:]
    lon[lon>179.99] = -360.+lon[lon>179.99]

    # Using the datetimes, split out the three time periods:
    # before, during, and after the event
    # ------------------------------------------------------
    print("Extracting times")
    total_times = np.array([base_date + relativedelta(days = ttime) for ttime in time])
    before_event = np.where((total_times >= dt_begin_start) & (total_times <= dt_begin_end))
    during_event = np.where((total_times >= dt_event_start) & (total_times <= dt_event_end))
    after_event  = np.where((total_times >= dt_end_start) & (total_times <= dt_end_end))

    # Create the HDF5 file and write the data
    # ---------------------------------------
    before_file = save_path + 'ceres_subset_' + case_date + '_before.hdf5'
    during_file = save_path + 'ceres_subset_' + case_date + '_during.hdf5'
    after_file  = save_path + 'ceres_subset_' + case_date + '_after.hdf5'

    before_dset = h5py.File(before_file,'w')
    during_dset = h5py.File(during_file,'w')
    after_dset  = h5py.File(after_file,'w')

    # Write the before data 
    print("Writing to out files")
    before_dset.create_dataset('time',        data = time[before_event])
    before_dset.create_dataset('latitude',    data = lat[before_event])
    before_dset.create_dataset('longitude',   data = lon[before_event])
    before_dset.create_dataset('swf',         data = swf[before_event])
    before_dset.create_dataset('lwf',         data = lwf[before_event])
    before_dset.create_dataset('sza',         data = sza[before_event])
    before_dset.create_dataset('alb',         data = alb[before_event])
    before_dset.create_dataset('clr_strong',  data = clr_str[before_event])
    before_dset.create_dataset('clr_weak',    data = clr_wek[before_event])

    # Write the during data 
    during_dset.create_dataset('time',        data = time[during_event])
    during_dset.create_dataset('latitude',    data = lat[during_event])
    during_dset.create_dataset('longitude',   data = lon[during_event])
    during_dset.create_dataset('swf',         data = swf[during_event])
    during_dset.create_dataset('lwf',         data = lwf[during_event])
    during_dset.create_dataset('sza',         data = sza[during_event])
    during_dset.create_dataset('alb',         data = alb[during_event])
    during_dset.create_dataset('clr_strong',  data = clr_str[during_event])
    during_dset.create_dataset('clr_weak',    data = clr_wek[during_event])

    # Write the after data 
    after_dset.create_dataset('time',        data = time[after_event])
    after_dset.create_dataset('latitude',    data = lat[after_event])
    after_dset.create_dataset('longitude',   data = lon[after_event])
    after_dset.create_dataset('swf',         data = swf[after_event])
    after_dset.create_dataset('lwf',         data = lwf[after_event])
    after_dset.create_dataset('sza',         data = sza[after_event])
    after_dset.create_dataset('alb',         data = alb[after_event])
    after_dset.create_dataset('clr_strong',  data = clr_str[after_event])
    after_dset.create_dataset('clr_weak',    data = clr_wek[after_event])

    data.close()
    
    # Save, write, and close the HDF5 file
    # --------------------------------------
    before_dset.close()
    during_dset.close()
    after_dset.close()

    print("Saved file", before_file)
    print("Saved file", during_file)
    print("Saved file", after_file)

# Data period is of format YYYYMMDDHH
def readgridCERES_hrly(data_dt,param,satellite = 'Aqua',minlat=60.0,\
        resolution = 0.25, season='all'):
    CERES_data_hrly = {}

    # Grab all the files
    if(satellite == 'Terra'):
        base_path = home_dir + '/data/CERES/SSF_Level2/Terra/'
    else:
        base_path = home_dir + '/data/CERES/SSF_Level2/Aqua/'
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
    
    # Set up values for gridding the flux data
    lat_gridder = minlat * (1. / resolution)
    max_idx = int(360 / resolution)
    multiplier = int(1 / resolution)
    adder = int(180 / resolution)
    
    lat_ranges = np.arange(minlat,90, resolution)
    lon_ranges = np.arange(-180,180, resolution)


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
            if(i % 10000 == 0):
                print(i)
            #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
            local_time = base_date + relativedelta(days = time[i])
            if((local_time >= dt_data_begin) & (local_time < dt_data_end)):
                if((flux[i] < 5000) and (flux[i] > 0) and (lat[i] >= minlat)):
                   #and (sza[i] < 80)):
                    # Skip over the really small flux values when doing
                    # daily averages
                    if(day and (sza[i] >= 77)):
                        continue 
                    ##!#index1 = int(np.floor(lat[i]*4 - lat_gridder))
                    ##!#index2 = int(np.floor(lon[i]*4 + 720.))
                    ##!#if(index1 < 0): index1 = 0
                    ##!#if(index1 > 719): index1 = 719
                    ##!#if(index2 < 0): index2 = 0                                                                                            
                    ##!#if(index2 > 1439): index2 = 1439

                    index1 = int(np.floor(lat[i]*multiplier - lat_gridder))
                    index2 = int(np.floor(lon[i]*multiplier + adder))
                    
                    if(index1 < 0): index1 = 0
                    if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                    if(index2 < 0): index2 = 0                                                                                            
                    if(index2 > (max_idx - 1)): index2 = max_idx - 1
             
                    try: 
                        swf_grid[index2, index1] = swf_grid[index2,index1] + flux[i]
                        #swf_grid[index2, index1] = (swf_grid[index2,index1]*count[index2,index1] + flux[i])/(count[index2,index1]+1)
                    except IndexError:
                        print(lat[i],lon[i],index1,index2)
                    #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                    #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                    #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                    count[index2, index1] = count[index2,index1] + 1
        data.close()
  
    swf_grid = swf_grid / count
 
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
    month_climo = np.zeros((6,CERES_data['data'].shape[1],CERES_data['data'].shape[2]))

    # Mask the monthly averages
    local_data = np.copy(CERES_data['data'][:,:,:])
    local_mask = np.ma.masked_where(local_data == -999., local_data)

    # Calculate monthly climatologies
    for m_i in range(6):
        month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::6,:,:],axis=0)
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

def plotCERES_hrly(pax, CERES_data_hrly, param, minlat=65, \
        vmin = None, vmax = None, title = '', label = None, \
        labelsize = None, labelticksize = 11, circle_bound = False, \
        gridlines = True, grid_data = True, zoom = True):

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
    if(label == None):
        label = param.upper() + ' [W/m$^{2}$]'
    pax.set_title(title)

    #plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
    mesh = pax.pcolormesh(plot_lon, plot_lat,mask_flux,transform = datacrs,\
        cmap = colormap, vmin = vmin, vmax = vmax, shading = 'auto')
    cbar = plt.colorbar(mesh,ax = pax, \
        extend = 'both', label = label)
    #cbar = plt.colorbar(mesh,ax = pax, orientation='vertical',\
    #    extend = 'both', fraction = 0.046, pad = 0.04)
    #cbar.set_label(label,fontsize = labelsize, weight='bold')
    #cbar.ax.tick_params(labelsize=labelticksize)

    pax.coastlines(resolution = '50m')
    #pax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)

    if(circle_bound):
        pax.set_boundary(circle, transform=pax.transAxes)

# Plot a daily average of CERES data
# ----------------------------------
def plotCERES_daily(CERES_data, pvar, end_str = None, satellite = 'Aqua',  \
        only_sea_ice = False, minlat = 65., avg_data = True, \
        vmin = None, vmax = None, cmap = 'viridis',\
        lat_circles = None, ax = None, save = False, min_ice = 0., \
        circle_bound = True, colorbar = True):

    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(isinstance(CERES_data, str)):
        date_str = CERES_data
        dt_date_str = datetime.strptime(date_str, '%Y%m%d')
        CERES_data = readgridCERES_daily(date_str, end_str = end_str, \
            satellite = satellite, minlat = minlat)
    else:
        dt_date_str = CERES_data['dt_dates'][0] 
        date_str = dt_date_str.strftime('%Y%m%d')

    # ----------------------------------------------------
    # Set up the overall figure
    # ----------------------------------------------------
    in_ax = True 
    if(ax is None): 
        in_ax = False
        plt.close('all')
        fig1 = plt.figure()
        #fig1 = plt.figure(figsize = (6,6))
        mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.Robinson()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    if(len(CERES_data[pvar].shape) == 3):
        pdata = np.nanmean(CERES_data[pvar], axis = 0)
        ice_data = np.nanmean(CERES_data['ice_conc'], axis = 0)
    else:
        pdata = CERES_data[pvar] 
        #ice_data = CERES_data['ice_conc']
        ice_data = np.nanmean(CERES_data['ice_conc'], axis = 0)

    pdata = np.ma.masked_where(ice_data < min_ice, pdata)

    mesh = ax.pcolormesh(CERES_data['lon'], CERES_data['lat'], pdata, \
        transform = ccrs.PlateCarree(), cmap = cmap, vmin = vmin, \
        vmax = vmax, shading = 'auto')
    ax.coastlines()
    ax.set_extent([-180, 180, minlat, 90], datacrs)
    
    if(colorbar):
        #cbar = plt.colorbar(mesh,ax = pax, orientation='horizontal',pad=0,\
        cbar = plt.colorbar(mesh,ax = ax, orientation='vertical',\
            pad = 0.04, fraction = 0.040)
        cbar.set_label(pvar, fontsize = 12, weight='bold')
        
    if(circle_bound):
        ax.set_boundary(circle, transform=ax.transAxes)

    plt.suptitle(date_str)

    # ----------------------------------------------------
    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax, lat_circles) 

    if(not in_ax):
        fig1.tight_layout()
        if(save):
            outname = 'ceres_' + satellite + '_daily_' + pvar + '_' + \
                date_str + '.png'
            fig1.savefig(outname, dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

def plotCERES_daily_allsat(date_str, pvar, end_str = None, \
        only_sea_ice = False, minlat = 65., avg_data = True, \
        lat_circles = None, save = False, \
        circle_bound = True, colorbar = True):

    # Read the data
    # -------------
    CERES_data1 = readgridCERES_daily(date_str, end_str = end_str, \
        satellite = 'Aqua',minlat = minlat)
    CERES_data2 = readgridCERES_daily(date_str, end_str = end_str, \
        satellite = 'Terra',minlat = minlat)
    CERES_data3 = readgridCERES_daily(date_str, end_str = end_str, \
        satellite = 'SuomiNPP',minlat = minlat)
    CERES_data4 = readgridCERES_daily(date_str, end_str = end_str, \
        satellite = 'SuomiNPP',minlat = minlat)

    total_data = np.concatenate((CERES_data1[pvar], CERES_data2[pvar], \
        CERES_data3[pvar]))
    total_data = np.ma.masked_where(total_data < 0, total_data)
    #total_data = np.nanmean(total_data, axis = 0)
    CERES_data4[pvar] = total_data 

    # Set up the overall figure
    # -------------------------
    fig = plt.figure(figsize = (9, 9))
    ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
    ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
    ax4 = fig.add_subplot(2,2,4, projection = mapcrs)

    # Plot the data
    # -------------
    plotCERES_daily(CERES_data1, pvar, end_str = end_str, satellite = 'Aqua',  \
        only_sea_ice = False, minlat = minlat, avg_data = True, \
        lat_circles = None, ax = ax1, save = False, \
        circle_bound = True, colorbar = True)
    plotCERES_daily(CERES_data2, pvar, end_str = end_str, satellite = 'Terra',  \
        only_sea_ice = False, minlat = minlat, avg_data = True, \
        lat_circles = None, ax = ax2, save = False, \
        circle_bound = True, colorbar = True)
    plotCERES_daily(CERES_data3, pvar, end_str = end_str, satellite = 'SuomiNPP',  \
        only_sea_ice = False, minlat = minlat, avg_data = True, \
        lat_circles = None, ax = ax3, save = False, \
        circle_bound = True, colorbar = True)
    plotCERES_daily(CERES_data4, pvar, end_str = end_str, satellite = 'SuomiNPP',  \
        only_sea_ice = False, minlat = minlat, avg_data = True, \
        lat_circles = None, ax = ax4, save = False, \
        circle_bound = True, colorbar = True)
    
    ax1.set_title('Aqua')
    ax2.set_title('Terra')
    ax3.set_title('SuomiNPP')
    ax4.set_title('Average')

    if(save):
        outname = 'ceres_allsat_' + pvar + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)

# Plot just a single swath on a 1-panel figure
#     skiprows - 
# NOTE: modified to work with readgridCERES_hrly_grid
# ---------------------------------------------------
def plotCERES_hrly_figure(date_str, param,  \
        only_sea_ice = False, minlat = 65., title = None, \
        lat_circles = None, grid_data = True, zoom = False, \
        vmin = None, vmax = None, circle_bound = True, \
        satellite = 'Aqua',
        ax = None, save = False):

    dt_date_str = datetime.strptime(date_str, "%Y%m%d%H%M")
    
    # ----------------------------------------------------
    # Read in data
    # ----------------------------------------------------
    if(grid_data):
        CERES_data_hrly = readgridCERES_hrly_grid(date_str[:10], param, \
            satellite = satellite, minlat = minlat)

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
    in_ax = True 
    if(ax is None): 
        in_ax = False
        plt.close('all')
        fig1 = plt.figure(figsize = (6,6))
        mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = init_proj(date_str)
        #mapcrs = ccrs.LambertConformal(central_longitude = -100.)
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    # ----------------------------------------------------
    # Use the single-swath plotting function to plot each
    # of the 3 data types
    # ----------------------------------------------------
    if(zoom):
        circle_bound = False
    plotCERES_hrly(ax, CERES_data_hrly, param, minlat = minlat, \
        vmin = vmin, vmax = vmax, title = title, label = 'Wm$^{-2}$', \
        circle_bound = circle_bound, gridlines = False, grid_data = grid_data)

    if(zoom):
        ax.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]['Lat'][1]],\
                        datacrs)
##!#        ax0.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
##!#                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
##!#                        datacrs)
    else:
        ax.set_extent([-180,180,minlat,90],ccrs.PlateCarree())

    #ax.set_title(date_str)

    # ----------------------------------------------------
    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax, lat_circles) 

    if(not in_ax):
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
    lon_ranges = CERES_data['lon'][0,:]
    #lon_ranges = np.arange(-180,180,1.0)

    # Make copy of CERES_data array
    print(CERES_data['dates'][month_idx::index_jumper])
    local_data   = np.copy(CERES_data['data'][month_idx::index_jumper,:,:])
    #local_counts = np.copy(CERES_data['OB_COUNT'][month_idx::index_jumper,:,:])
    #local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where(((local_data == -999.) | \
        (CERES_data['lat'] < minlat)), local_data)
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
                if(trend_type=='standard'): 
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

    #ceres_trends = np.ma.masked_where(((CERES_data['lat'] < minlat) | \
    #    (ceres_trends == -999.)), ceres_trends)
    ceres_trends = np.ma.masked_where(CERES_data['lat'] < minlat, ceres_trends)
    ceres_pvals  = np.ma.masked_where(CERES_data['lat'] < minlat, ceres_pvals)
    ceres_uncert = np.ma.masked_where(CERES_data['lat'] < minlat, ceres_uncert)

    return ceres_trends, ceres_pvals, ceres_uncert

# title is the plot title
# ptype is 'climo', 'trend'
def plotCERES_spatial(pax, plat, plon, pdata, ptype, ptitle = '', plabel = '', \
        vmin = None, vmax = None, colorbar_label_size = 16, minlat = 65., \
        colorbar = True, \
        pvals = None):

    if(vmin == None):
        vmin = np.nanmin(pdata)
    if(vmax == None):
        vmax = np.nanmax(pdata)

    if(ptype == 'trend'):
        colormap = plt.cm.bwr
    elif(ptype == 'uncert'):
        #colormap = plt.cm.plasma
        colormap = plt.cm.get_cmap('jet', 6)
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

    #pax.gridlines()
    pax.coastlines(resolution='50m')
    mesh = pax.pcolormesh(plon, plat,\
            mask_AI.T,transform = datacrs,shading = 'auto',\
            cmap = colormap,vmin=vmin,vmax=vmax)
    if(pvals is not None):
        cyclic_pvals, cyclic_lons = add_cyclic_point(pvals, plon[0,:])
        print(pvals.shape, plat.shape, plat2.shape)
        mask_pvals = np.ma.masked_where((plat < minlat) | \
            (pvals > 0.05), pvals)
        pax.pcolor(plon, plat, mask_pvals, hatch = '...', alpha = 0.0, \
            shading = 'auto', transform = datacrs)

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
    if(colorbar):
        cbar = plt.colorbar(mesh,ticks = np.arange(0.0,400.,25.),orientation='horizontal',pad=0,\
            aspect=50,shrink = 0.905,label=CERES_data['parm_name'])
        #cbar = plt.colorbar(mesh,\
        #    ax = pax, orientation='vertical',shrink = 0.8, extend = 'both')
        ##cbar10.set_label('UV Aerosol Index',weight='bold',fontsize=colorbar_label_size)
        ##cbar.ax.tick_params(labelsize=14)
        #cbar.set_label(plabel,fontsize=colorbar_label_size,weight='bold')
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
def plotCERES_MonthClimo(CERES_data,month_idx,minlat = 60, ax = None, \
        colorbar = True, title = None):

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
    if(title == None):
        title = 'CERES ' + CERES_data['param'] + ' ' + str_month + \
            ' Climatology\n'+str_month+'. ' + CERES_data['dates'][0][:4] \
            + ' - '+str_month+'. ' + CERES_data['dates'][-1][:4]

    # Make figure
    in_ax = True 
    if(ax is None): 
        in_ax = False
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
    if(colorbar):
        #cbar = plt.colorbar(mesh,ax = pax, orientation='horizontal',pad=0,\
        cbar = plt.colorbar(mesh,ax = ax, orientation='vertical',\
            pad = 0.04, fraction = 0.040)
        cbar.set_label('W/m2', fontsize = 12, weight='bold')
    #cbar = plt.colorbar(mesh,ticks = np.arange(0.0,400.,25.),orientation='horizontal',pad=0,\
    #    aspect=50,shrink = 0.845,label=CERES_data['parm_name'])
    ax.set_title(title)

    if(not in_ax):
        fig1.tight_layout()
        if(save):
            outname = 'ceres_' + CERES_data['satellite'] + '_daily_' + pvar + '_' + \
                date_str + '_FINISH.png'
            fig1.savefig(outname, dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

# Designed to work with the netCDF data
def plotCERES_MonthTrend(CERES_data,month_idx=None,save=False,\
        trend_type='standard',season='',minlat=65.,return_trend=False, \
        colorbar = True, colorbar_label_size = None,title = None, \
        pax = None, show_pval = False, uncert_ax = None):

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
    ceres_trends, ceres_pvals, ceres_uncert = \
        calcCERES_grid_trend(CERES_data, month_idx, trend_type, \
        minlat)

    if(not show_pval):
        ceres_pvals = None
    else:
        print('month_idx = ',month_idx,' PVAL nanmean = ', \
            np.nanmean(ceres_pvals))

    if(uncert_ax is None):
        ceres_uncert = None
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

    if(title is None):
        title = 'CERES ' + CERES_data['param'] + ' ' + month_string + 'Trends'\
            '\n'+start_date.strftime('%b. %Y') + ' - ' +\
            end_date.strftime('%b. %Y')

    # Call plotCERES_spatial to add the data to the figure

    ii = 15
    jj = 180
    print(CERES_data['param'])
    print(CERES_data['lat'][ii,jj], CERES_data['lon'][ii,jj])
    print(CERES_data['data'][month_idx::index_jumper,ii,jj])

    if(pax is None):
        plt.close('all')
        fig1 = plt.figure(figsize = (6,6))
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

        plotCERES_spatial(ax, CERES_data['lat'], CERES_data['lon'], \
            ceres_trends, 'trend', ptitle = title, plabel = 'W m$^{-2}$ per study period', \
            vmin = v_min, vmax = v_max, colorbar_label_size = colorbar_label_size, \
            minlat = minlat, pvals = ceres_pvals)

        fig1.tight_layout()

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
            ceres_trends, 'trend', ptitle = title, plabel = 'W m$^{-2}$ per study period', \
            vmin = v_min, vmax = v_max, colorbar_label_size = colorbar_label_size, \
            minlat = minlat)

    if(uncert_ax is not None):
        plotCERES_spatial(uncert_ax, CERES_data['lat'], CERES_data['lon'], \
            ceres_uncert, 'uncert', ptitle = title, plabel = 'W/m2', \
            colorbar = colorbar, colorbar_label_size = colorbar_label_size, \
            vmin = 0, vmax = 20.0, minlat = minlat)

    if(return_trend == True):
        return ceres_trends

def plotCERES_MonthTrend_AllClr(start_date, end_date, month_idx, \
        minlat = 65., satellite = 'Aqua', save = False):

    CERES_swall = readgridCERES(start_date,end_date,'toa_sw_all_mon',\
        satellite = satellite,minlat=minlat,calc_month = True,season = 'sunlight')
    CERES_lwall = readgridCERES(start_date,end_date,'toa_lw_all_mon',\
        satellite = satellite,minlat=minlat,calc_month = True,season = 'sunlight')
    CERES_swclr = readgridCERES(start_date,end_date,'toa_sw_clr_mon',\
        satellite = satellite,minlat=minlat,calc_month = True,season = 'sunlight')
    CERES_lwclr = readgridCERES(start_date,end_date,'toa_lw_clr_mon',\
        satellite = satellite,minlat=minlat,calc_month = True,season = 'sunlight')
    #CERES_all  = readgridCERES(start_date,end_date,'toa_sw_all_mon',\
    #    satellite = 'Aqua',minlat=60.5,calc_month = True,season = 'sunlight')
    
    fig1 = plt.figure(figsize = (10,10))
    ax1 = fig1.add_subplot(2,2,1, projection = mapcrs)
    ax2 = fig1.add_subplot(2,2,2, projection = mapcrs)
    ax3 = fig1.add_subplot(2,2,3, projection = mapcrs)
    ax4 = fig1.add_subplot(2,2,4, projection = mapcrs)
    
    plotCERES_MonthTrend(CERES_swclr,month_idx=month_idx,save=False,\
        trend_type='standard',season='sunlight',minlat=minlat,return_trend=False, \
        pax = ax1)
    plotCERES_MonthTrend(CERES_lwclr,month_idx=month_idx,save=False,\
        trend_type='standard',season='sunlight',minlat=minlat,return_trend=False, \
        pax = ax2)
    plotCERES_MonthTrend(CERES_swall,month_idx=month_idx,save=False,\
        trend_type='standard',season='sunlight',minlat=minlat,return_trend=False, \
        pax = ax3)
    plotCERES_MonthTrend(CERES_lwall,month_idx=month_idx,save=False,\
        trend_type='standard',season='sunlight',minlat=minlat,return_trend=False, \
        pax = ax4)
    
    fig1.tight_layout()

    if(save):
        outname = '_'.join(['ceres',satellite,'month'+str(month_idx),\
            'AllClr']) + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)

    else:
        plt.show()

# Generate a 15-panel figure comparing the climatology and trend between 3
# versions of the CERES data for all months
def plotCERES_ClimoTrend_all(CERES_data,\
        trend_type = 'standard', minlat=65.,save=False):

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
    fig = plt.figure(figsize=(9.2,13))
    #plt.suptitle('NAAPS Comparisons: '+start_date.strftime("%B"),y=0.95,\
    #    fontsize=18,fontweight=4,weight='bold')
    gs = gridspec.GridSpec(nrows=6, ncols=3, hspace = 0.001, wspace = 0.15)

    # - - - - - - - - - - - - - - - - - - - - -
    # Plot the climatologies along the top row
    # - - - - - - - - - - - - - - - - - - - - -
       
    # Plot DATA1 climos
    # -----------------
    ##!## Make copy of CERES_data array
    local_data1_Apr  = np.copy(CERES_data['MONTH_CLIMO'][0,:,:])
    local_data1_May  = np.copy(CERES_data['MONTH_CLIMO'][1,:,:])
    local_data1_Jun  = np.copy(CERES_data['MONTH_CLIMO'][2,:,:])
    local_data1_Jul  = np.copy(CERES_data['MONTH_CLIMO'][3,:,:])
    local_data1_Aug  = np.copy(CERES_data['MONTH_CLIMO'][4,:,:])
    local_data1_Sep  = np.copy(CERES_data['MONTH_CLIMO'][5,:,:])

    mask_AI1_Apr = np.ma.masked_where(local_data1_Apr == -999.9, local_data1_Apr)
    mask_AI1_May = np.ma.masked_where(local_data1_May == -999.9, local_data1_May)
    mask_AI1_Jun = np.ma.masked_where(local_data1_Jun == -999.9, local_data1_Jun)
    mask_AI1_Jul = np.ma.masked_where(local_data1_Jul == -999.9, local_data1_Jul)
    mask_AI1_Aug = np.ma.masked_where(local_data1_Aug == -999.9, local_data1_Aug)
    mask_AI1_Sep = np.ma.masked_where(local_data1_Sep == -999.9, local_data1_Sep)
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,CERES_data1['LON'][0,:])
    ##!#plat,plon = np.meshgrid(CERES_data1['LAT'][:,0],cyclic_lons)   
  
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
    plotCERES_MonthClimo(CERES_data,0,minlat = minlat, ax = ax00, title = '')
    plotCERES_MonthTrend(CERES_data,month_idx=0,save=False,\
        trend_type='standard',season='sunlight',minlat=minlat,\
        return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax01, show_pval = True, uncert_ax = ax02, title = '')
    #plotCERES_spatial(ax00, CERES_data['lats'], CERES_data['lons'], mask_AI1_Apr, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 1.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotCERES_MonthTrend(CERES_data,month_idx=0,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax01, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax02)
    #plotCERES_MonthTrend(CERES_data,month_idx=month_idx,save=False,\
    #    trend_type='standard',season='sunlight',minlat=65.,return_trend=False, \
    #    pax = ax1)

    # Plot the figures in the first row: May
    # ---------------------------------------
    #plotCERES_spatial(ax10, CERES_data['lats'], CERES_data['lons'], mask_AI1_May, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 2.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotCERES_MonthTrend(CERES_data,month_idx=1,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax11, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax12)
    plotCERES_MonthClimo(CERES_data,1,minlat = minlat, ax = ax10, title = '')
    plotCERES_MonthTrend(CERES_data,month_idx=1,save=False,\
        trend_type='standard',season='sunlight',\
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax11, show_pval = True, uncert_ax = ax12, title = '')

    # Plot the figures in the first row: June
    # ---------------------------------------
    #plotCERES_spatial(ax20, CERES_data['lats'], CERES_data['lons'], mask_AI1_Jun, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotCERES_MonthTrend(CERES_data,month_idx=2,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax21, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax22)
    plotCERES_MonthClimo(CERES_data,2,minlat = minlat, ax = ax20, title = '')
    plotCERES_MonthTrend(CERES_data,month_idx=2,save=False,\
        trend_type='standard',season='sunlight',
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax21, show_pval = True, uncert_ax = ax22, title = '')

    # Plot the figures in the second row: July
    # ----------------------------------------
    #plotCERES_spatial(ax30, CERES_data['lats'], CERES_data['lons'], mask_AI1_Jul, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotCERES_MonthTrend(CERES_data,month_idx=3,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax31, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax32)
    plotCERES_MonthClimo(CERES_data,3,minlat = minlat, ax = ax30, title = '')
    plotCERES_MonthTrend(CERES_data,month_idx=3,save=False,\
        trend_type='standard',season='sunlight',\
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax31, show_pval = True, uncert_ax = ax32, title = '')

    # Plot the figures in the third row: August
    # -----------------------------------------
    #plotCERES_spatial(ax40, CERES_data['lats'], CERES_data['lons'], mask_AI1_Aug, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotCERES_MonthTrend(CERES_data,month_idx=4,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax41, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax42)
    plotCERES_MonthClimo(CERES_data,4,minlat = minlat, ax = ax40, title = '')
    plotCERES_MonthTrend(CERES_data,month_idx=4,save=False,\
        trend_type='standard',season='sunlight',\
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax41, show_pval = True, uncert_ax = ax42, title = '')

    # Plot the figures in the third row: September
    # --------------------------------------------
    #plotCERES_spatial(ax50, CERES_data['lats'], CERES_data['lons'], mask_AI1_Sep, \
    #    'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 1.0, \
    #    minlat = minlat, colorbar = cbar_switch)
    #plotCERES_MonthTrend(CERES_data,month_idx=5,trend_type=trend_type,label = ' ',\
    #    minlat=65.,title = ' ', pax = ax51, colorbar = cbar_switch, \
    #    colorbar_label_size = colorbar_label_size, show_pval = True, \
    #    uncert_ax = ax52)
    plotCERES_MonthClimo(CERES_data,5,minlat = minlat, ax = ax50, title = '')
    plotCERES_MonthTrend(CERES_data,month_idx=5,save=False,\
        trend_type='standard',season='sunlight',\
        minlat=minlat,return_trend=False, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax51, show_pval = True, uncert_ax = ax52, title = '')

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

    plt.suptitle(CERES_data['param'])

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

    outname = '_'.join(['ceres',CERES_data['param'],'grid','comps','all6']) + '.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()


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

##!#def plot_compare_OMI_CERES_hrly(OMI_date,CERES_date,minlat=65,max_AI = -200.,\
##!#        omi_dtype = 'control', only_sea_ice = False, skiprows = None, save=False):
##!##def plot_compare_OMI_CERES_hrly(OMI_hrly,CERES_hrly,minlat=65,save=False):
##!#
##!#    if('/home/bsorenson/Research/OMI' not in sys.path):
##!#        sys.path.append('/home/bsorenson/Research/OMI')
##!#    from OMILib import readOMI_swath_shawn, readOMI_swath_hdf, \
##!#        plotOMI_single_swath
##!#
##!#    # ----------------------------------------------------
##!#    # Set up the overall figure
##!#    # ----------------------------------------------------
##!#    plt.close('all')
##!#    fig1 = plt.figure(figsize = (6,6))
##!#    mapcrs = ccrs.NorthPolarStereo()
##!#    ax0 = fig1.add_subplot(1,1,1, projection = mapcrs)
##!#
##!#    # Set up mapping variables 
##!#    datacrs = ccrs.PlateCarree() 
##!#    colormap = plt.cm.jet
##!#    if(minlat < 45):
##!#        mapcrs = ccrs.Miller()
##!#    else:
##!#        mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)
##!#
##!#    # Step 1: read in OMI and CERES data
##!#    # ----------------------------------
##!#    #OMI_hrly = readOMI_single_swath(OMI_date,60,only_sea_ice = only_sea_ice)
##!#    if(omi_dtype == 'shawn'):
##!#        OMI_hrly  = readOMI_swath_shawn(OMI_date, latmin = minlat,\
##!#            skiprows = skiprows)
##!#    else:
##!#        if(omi_dtype == 'jz'):
##!#            omi_dtype = 'JZ'
##!#        OMI_hrly  = readOMI_swath_hdf(OMI_date, omi_dtype, \
##!#            only_sea_ice = only_sea_ice, latmin = minlat, \
##!#            skiprows = skiprows)
##!#
##!#    CERES_hrly = readgridCERES_hrly(CERES_date,'SWF')
##!#
##!#    # Step 2: Set up figure to have 3 panels
##!#    # --------------------------------------
##!#    plt.close('all')
##!#    #fig, axs = plt.subplots(1,3,subplot_kw={'projection': mapcrs})
##!#    fig = plt.figure(figsize=(17,5))
##!#    ax0 = fig.add_subplot(1,3,2,projection = mapcrs)
##!#    ax1 = fig.add_subplot(1,3,3,projection = mapcrs)
##!#    ax2 = fig.add_subplot(1,3,1)
##!#
##!#    # Step 3: Plot OMI and CERES data in first two panels
##!#    # ---------------------------------------------------
##!#
##!#    # -------------------------------------------------------
##!#    # Use the single-swath plotting function to plot OMI data
##!#    # -------------------------------------------------------
##!#    plotOMI_single_swath(ax0, OMI_hrly, title = omi_dtype.title(), \
##!#        circle_bound = True, gridlines = False)
##!#
##!#    ax0.set_extent([-180,180,minlat,90], datacrs)
##!#
##!#
##!#    local_data = np.copy(CERES_hrly['data'])
##!#
##!#    # -------------------------------------------------------
##!#    # Use the single-swath plotting function to plot OMI data
##!#    # -------------------------------------------------------
##!#    plot_lat, plot_lon = np.meshgrid(CERES_hrly['lat'],CERES_hrly['lon'])
##!#    mask_flux = np.ma.masked_where(((CERES_hrly['counts'] == 0) | \
##!#        (CERES_hrly['lat'] < minlat)), local_data)
##!#
##!#    ax1.gridlines(ylocs = np.arange(minlat,90,5))
##!#    ax1.coastlines(resolution = '50m')
##!#    ax1.set_title('CERES ' + CERES_hrly['param'] + ' ' + CERES_hrly['date'])
##!#    #plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
##!#    mesh1 = ax1.pcolormesh(plot_lon, plot_lat,mask_flux,transform = datacrs,cmap = colormap,\
##!#            vmin = min_dict[CERES_hrly['param']], \
##!#            vmax = max_dict[CERES_hrly['param']])
##!#    ax1.set_extent([-180,180,minlat,90],datacrs)
##!#    ax1.set_boundary(circle, transform=ax1.transAxes)
##!#            #vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])
##!#    cbar = plt.colorbar(mesh1,ax=ax1,orientation='vertical',pad=0.02,\
##!#        aspect=50,shrink = 1.000,label=CERES_hrly['param'] + ' [W/m2]')
##!#    
##!#    # Step 4: Plot scatter OMI and CERES comparison in third panel
##!#    # ------------------------------------------------------------
##!#
##!#    # Make gridded lat and lon arrays to use for coloring the plot
##!#    glat, glon = np.meshgrid(OMI_hrly['LAT'],OMI_hrly['LON'])
##!#    
##!#    # Mask any empty boxes
##!#    mask_AI = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
##!#        (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),OMI_hrly['AI'])
##!#    mask_flux = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
##!#        (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),CERES_hrly['data'])
##!#    mask_sza  = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
##!#        (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),OMI_hrly['SZA'])
##!#    #mask_lat  = np.ma.masked_where(((OMI_hrly['AI_count'] == 0) | \
##!#    #    (CERES_hrly['counts'] == 0) | (OMI_hrly['AI'] < max_AI)),glat)
##!#    ##!## Convert the index to a string using datetime
##!#    ##!#if(month != None):
##!#    ##!#    dt_obj = datetime.strptime(OMI_data['DATES'][month],"%Y%m")
##!#    ##!#    title = 'OMI AI / CERES ' + CERES_data['param'] + '\n'+ \
##!#    ##!#        dt_obj.strftime("%b") + " Trend Comparison"
##!#    ##!#    outname = 'omi_ceres_trend_comp_'+dt_obj.strftime("%b")+'_'+\
##!#    ##!#        OMI_data['VERSION']+'vCERES.png'
##!#    ##!#else:
##!#    ##!#    title = 'OMI AI / CERES ' + CERES_data['param'] + ' Trend Comparison'
##!#    ##!#    outname = 'omi_ceres_trend_comp_'+\
##!#    ##!#        OMI_data1['VERSION']+'vCERES.png'
##!#
##!#
##!#    #print("Pearson:  ",pearsonr(mask_AI,mask_flux))
##!#    #print("Spearman: ",spearmanr(mask_AI,mask_flux))
##!#
##!#    ##!#xy = np.vstack([mask_AI,mask_flux])
##!#    ##!#z = stats.gaussian_kde(xy)(xy)
##!#
##!#    ##!## Plot a somewhat-robust best fit line using Huber Regression
##!#    ##!## -----------------------------------------------------------
##!#    ##!#x_scaler,y_scaler = StandardScaler(), StandardScaler()
##!#    ##!#x_train = x_scaler.fit_transform(mask_AI[...,None])
##!#    ##!#y_train = y_scaler.fit_transform(mask_flux[...,None])
##!#
##!#    ##!#model = HuberRegressor(epsilon=1)
##!#    ##!#model.fit(x_train,y_train)
##!#    ##!#
##!#    ##!#test_x = np.array([np.min(mask_AI),np.max(mask_AI)])
##!#    ##!#predictions = y_scaler.inverse_transform(\
##!#    ##!#    model.predict(x_scaler.transform(test_x[...,None])))
##!#
##!#    ##!## One to one line stuff
##!#    ##!#xs = np.arange(np.min(mask_AI),np.max(mask_AI),0.1)
##!#
##!#    #plt.close('all')
##!#    #fig1 = plt.figure()
##!#    #scat = ax2.scatter(mask_AI,mask_flux,c=mask_lat,s=8)
##!#    scat = ax2.scatter(mask_AI,mask_flux,c=mask_sza,s=8)
##!#    cbar = plt.colorbar(scat,ax=ax2,orientation='vertical',\
##!#        label='Solar Zenith Angle',pad=0.02,aspect=50)
##!#    ax2.set_title(OMI_hrly['date'])
##!#    #plt.title(OMI_hrly['date'])
##!#    #plt.scatter(mask_AI,mask_flux,c=z,s=8)
##!#    #plt.plot(test_x,predictions,color='tab:green',linestyle='--',label='Huber Fit')
##!#    #plt.plot(test_x,predictions,color='tab:green',linestyle='--',label='Huber Fit')
##!#    # Plot an unrobust fit line using linear regression
##!#    # -------------------------------------------------
##!#    #plt.plot(np.unique(mask_AI),np.poly1d(np.polyfit(mask_AI,\
##!#    #    mask_flux,1))(np.unique(mask_AI)),color='tab:orange',\
##!#    #    linestyle='--',label='Polyfit Fit')
##!#    # Plot a one-to-one line
##!#    #plt.plot(xs,xs,label='1-1',color='tab:red')
##!#    ##if((month == 0) | (month == 1)):
##!#    ##    plt.xlim(-0.5,0.3)
##!#    ##    plt.ylim(-0.5,0.3)
##!#    ##elif((month == 2)):
##!#    ##    plt.xlim(-0.6,0.5)
##!#    ##    plt.ylim(-0.6,0.5)
##!#    ##elif((month == 3)):
##!#    ##    plt.xlim(-0.5,0.5)
##!#    ##    plt.ylim(-0.5,0.5)
##!#    ##elif((month == 4)):
##!#    ##    plt.xlim(-0.5,0.7)
##!#    ##    plt.ylim(-0.5,0.7)
##!#    ##elif((month == 5)):
##!#    ##    plt.xlim(-0.5,0.5)
##!#    ##    plt.ylim(-0.5,0.5)
##!#    ##else:
##!#    ##    plt.xlim(-0.3,0.3)
##!#    ##    plt.ylim(-0.3,0.3)
##!#    #axs[2].legend()
##!#    ax2.set_xlabel('OMI AI')
##!#    ax2.set_ylabel('CERES ' + CERES_hrly['param'])
##!#    #plt.subplots_adjust(wspace=0.1, hspace=0.10)
##!#    plt.tight_layout()
##!#    #plt.title(title)
##!#    outname = 'omi_ceres_compare_'+OMI_date+'.png'
##!#    if(save == True):
##!#        plt.savefig(outname,dpi=300)
##!#        print("Saved image",outname)
##!#    else:
##!#        plt.show()
##!#    return mask_AI,mask_flux

#def plot_compare_OMI_CERES_hrly_grid(OMI_date,CERES_date,minlat=65,max_AI = -200.,\
def plot_compare_OMI_CERES_hrly_grid(date_str,minlat=65,max_AI = -200.,\
        omi_dtype = 'control', ceres_dtype = 'swf', only_sea_ice = False, \
        only_ice = False, no_ice = False,  skiprows = None, save=False, \
        composite = False, zoom = False, lat_circles = None, \
        show_scatter = False):
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
    print(CERES_hrly.keys())
    if(CERES_hrly == -1):
        print("ERROR: invalid value returned from readgridCERES_hrly_grid")
        print("returning")
        return

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
    #mapcrs = init_proj(date_str)
    if(show_scatter):
        fig = plt.figure(figsize=(9,6))
        ax0 = fig.add_subplot(2,2,2,projection = mapcrs)
        ax1 = fig.add_subplot(2,2,3,projection = mapcrs)
        ax3 = fig.add_subplot(2,2,1)
        #ax3 = fig.add_subplot(2,2,4,projection = crs1)
        ax2 = fig.add_subplot(2,2,4)
    else:
        fig = plt.figure(figsize=(14,3))
        ax0 = fig.add_subplot(1,3,2,projection = mapcrs)
        ax1 = fig.add_subplot(1,3,3,projection = mapcrs)
        ax3 = fig.add_subplot(1,3,1)
        #ax3 = fig.add_subplot(2,2,4,projection = crs1)

    ##!## Plot the true-color data for the previous date
    ##!## ----------------------------------------------
    ##!#ax3.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
    ##!#    var1.y[-1], var1.y[0]), origin='upper')

    ##!## Zoom in the figure if desired
    ##!## -----------------------------
    ##!#if(zoom):
    ##!#    ax3.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
    ##!#        lat_lims1[1]],crs = datacrs)

    # Read in the matching Worldview image
    # ------------------------------------
    base_path = '/home/bsorenson/Research/OMI/'
    img_name = dt_date_str.strftime(base_path + 'snapshot-%Y-%m-%dT00_00_00Z_zoomed2.jpg') 
    print("Looking for ",img_name)

    pic = plt.imread(img_name)
    ax3.imshow(pic)
    ax3.axis('off')

    # ------------------------------------------------------------------------
    #
    # Step 3: Plot OMI and CERES data in first two panels
    #
    # ------------------------------------------------------------------------
    lat_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                dt_date_str.strftime('%H%M')]['Lat']
    lon_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                dt_date_str.strftime('%H%M')]['Lon']
    ##!#if(zoom):
    ##!#    OMI_base['UVAI'] = np.ma.masked_where(~((OMI_base['LAT'] >= lat_lims[0]) & \
    ##!#                                            (OMI_base['LAT'] <  lat_lims[1]) & \
    ##!#                                            (OMI_base['LON'] >= lon_lims[0]) & \
    ##!#                                            (OMI_base['LON'] <  lon_lims[1])), OMI_base['UVAI'])
    ##!#    #OMI_base['UVAI'] = np.ma.masked_where(OMI_base['UVAI'] > 1, OMI_base['UVAI'])
    ##!#    CERES_hrly[ceres_dtype] = np.ma.masked_where(~((CERES_hrly['lat'] >= lat_lims[0]) & \
    ##!#                                              (CERES_hrly['lat'] <  lat_lims[1]) & \
    ##!#                                              (CERES_hrly['lon'] >= lon_lims[0]) & \
    ##!#                                              (CERES_hrly['lon'] <  lon_lims[1])), CERES_hrly[ceres_dtype])

    ##!## Extract OMI and CERES values at a desired point
    ##!## -----------------------------------------------
    ##!#tlat = 78.4375
    ##!#tlon = 138.9517
    ##!##tlat = 85.1155
    ##!##tlon = -68.2791
    ##!#c_idx = nearest_gridpoint(tlat, tlon,\
    ##!#    CERES_hrly['lat'], CERES_hrly['lon'])
    ##!#if(len(c_idx[0]) > 1):
    ##!#    c_idx = (np.array([c_idx[0][0]])), (np.array([c_idx[1][0]]))
    ##!#tflux = CERES_hrly[ceres_dtype][c_idx]
    ##!#tfluxlat = CERES_hrly['lat'][c_idx]
    ##!#tfluxlon = CERES_hrly['lon'][c_idx]

    ##!#o_idx = nearest_gridpoint(tlat, tlon,\
    ##!#    OMI_base['LAT'], OMI_base['LON'])
    ##!#if(len(o_idx[0]) > 1):
    ##!#    print('o_idx = ',o_idx)
    ##!#    o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
    ##!#tai    = OMI_base['UVAI'][o_idx]
    ##!#tailat = OMI_base['LAT'][o_idx]
    ##!#tailon = OMI_base['LON'][o_idx]

    ##!#ax0.plot(tailon, tailat,
    ##!#         color='k', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax0.plot(tfluxlon, tfluxlat,
    ##!#         color='r', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax0.plot(tlon, tlat,
    ##!#         color='w', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax1.plot(tailon, tailat,
    ##!#         color='k', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax1.plot(tfluxlon, tfluxlat,
    ##!#         color='r', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax1.plot(tlon, tlat,
    ##!#         color='w', linewidth=2, marker='o',
    ##!#         transform=datacrs)

    ##!#print("Selected OMI and CERES at ", tlat, ' ', tlon,':')
    ##!#print(' OMI   - ',tai)
    ##!#print('     lat - ', tailat)
    ##!#print('     lon - ', tailon)
    ##!#print(' CERES - ',tflux)
    ##!#print('     lat - ', tfluxlat)
    ##!#print('     lon - ', tfluxlon)

    # Find averages of the data within the ergion


    # Use the single-swath plotting function to plot OMI data
    # -------------------------------------------------------
    plotOMI_single_swath(ax0, OMI_base, title = None, \
        circle_bound = circle_bound, gridlines = False)
    #plotOMI_single_swath(ax0, OMI_base, title = 'OMI ' + omi_dtype + \
    #    ' ' + OMI_date, circle_bound = circle_bound, gridlines = False)

    # Use the single-swath plotting function to plot CERES data
    # ---------------------------------------------------------
    plotCERES_hrly(ax1, CERES_hrly, 'swf', minlat = minlat, \
        vmin = 90, vmax = 250, title = '', label = 'TOA Flux [Wm$^{-2}$]', \
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

    labelsize = 13
    scat_add = ''
    if(show_scatter):
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
            ceres_match_flux[ii] = CERES_hrly[ceres_dtype][o_idx]
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

        scat_add = '_scatter'
        # Mask any empty boxes
        print('Before scatter', np.ma.masked_invalid(ceres_match_flux).compressed().shape)
        scat = ax2.scatter(mask_UVAI,ceres_match_flux,s=8)
        #scat = ax2.scatter(mask_AI,mask_flux,c=mask_sza,s=8)
        #cbar = plt.colorbar(scat,ax=ax2,orientation='vertical',\
        #    label='Viewing Zenith Angle',pad=0.02,aspect=50)
        #ax2.set_title(OMI_base['date'])
        ax2.set_xlabel('OMI AI', weight = 'bold')
        ax2.set_ylabel('CERES ' + CERES_hrly['param'], weight = 'bold')
        plot_subplot_label(ax2, '(d)', fontsize = labelsize, backgroundcolor = 'white')

    fig.tight_layout()

    plot_figure_text(ax3, 'MODIS True Color', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax0, 'OMI AI', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax1, 'CERES SWF', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax3, '(a)', fontsize = labelsize, backgroundcolor = 'white')
    plot_subplot_label(ax0, '(b)', fontsize = labelsize, backgroundcolor = 'white')
    plot_subplot_label(ax1, '(c)', fontsize = labelsize, backgroundcolor = 'white')

    if(save == True):
        outname = 'omi_ceres_compare_'+OMI_date+scat_add+'.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

    #return OMI_base, CERES_hrly

def plot_compare_OMI_CERES_hrly_grid_2case(date_str1,date_str2, minlat=65,max_AI = -200.,\
        omi_dtype = 'control', ceres_dtype = 'swf', only_sea_ice = False, \
        only_ice = False, no_ice = False,  skiprows = None, save=False, \
        composite = False, zoom = False, lat_circles = None, \
        show_scatter = False):
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
    dt_date_str1 = datetime.strptime(date_str1,"%Y%m%d%H%M")
    dt_date_str2 = datetime.strptime(date_str2,"%Y%m%d%H%M")
    OMI_date1 = aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][\
                   dt_date_str1.strftime('%H%M')]['omi']
    OMI_date2 = aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][\
                   dt_date_str2.strftime('%H%M')]['omi']
    OMI_date1 = OMI_date1.split('_')[-2][:4] + OMI_date1.split('_')[-2][5:9] + \
        OMI_date1.split('_')[-2][10:14]
    OMI_date2 = OMI_date2.split('_')[-2][:4] + OMI_date2.split('_')[-2][5:9] + \
        OMI_date2.split('_')[-2][10:14]
    CERES_date1 = dt_date_str1.strftime('%Y%m%d') + \
        aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][\
        dt_date_str1.strftime('%H%M')]['ceres_time']
    CERES_date2 = dt_date_str2.strftime('%Y%m%d') + \
        aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][\
        dt_date_str2.strftime('%H%M')]['ceres_time']

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
        OMI_base1  = readOMI_swath_shawn(OMI_date1, latmin = minlat)
        OMI_base2  = readOMI_swath_shawn(OMI_date2, latmin = minlat)
    else:
        OMI_base1  = readOMI_swath_hdf(OMI_date1, omi_dtype, \
            only_sea_ice = only_sea_ice, only_ice = only_ice, \
            no_ice = no_ice, latmin = minlat, \
            skiprows = skiprows)
        OMI_base2  = readOMI_swath_hdf(OMI_date2, omi_dtype, \
            only_sea_ice = only_sea_ice, only_ice = only_ice, \
            no_ice = no_ice, latmin = minlat, \
            skiprows = skiprows)

    CERES_hrly1 = readgridCERES_hrly_grid(CERES_date1, 'SWF', minlat = 55.)
    if(CERES_hrly1 == -1):
        print("ERROR: invalid value returned from readgridCERES_hrly_grid")
        print("returning")
        return
    CERES_hrly2 = readgridCERES_hrly_grid(CERES_date2, 'SWF', minlat = 55.)
    if(CERES_hrly2 == -1):
        print("ERROR: invalid value returned from readgridCERES_hrly_grid")
        print("returning")
        return

    # ------------------------------------------------------------------------
    #
    # Step 2: Set up figure to have 6 panels
    #
    # ------------------------------------------------------------------------
    plt.close('all')
    #fig, axs = plt.subplots(1,3,subplot_kw={'projection': mapcrs})
    mapcrs1 = init_proj(date_str1)
    mapcrs2 = ccrs.NorthPolarStereo(central_longitude = -40)
    fig = plt.figure(figsize=(14,6))
    ax0 = fig.add_subplot(2,3,1)
    ax1 = fig.add_subplot(2,3,2,projection = mapcrs1)
    ax2 = fig.add_subplot(2,3,3,projection = mapcrs1)
    ax3 = fig.add_subplot(2,3,4)
    ax4 = fig.add_subplot(2,3,5,projection = mapcrs2)
    ax5 = fig.add_subplot(2,3,6,projection = mapcrs2)
    #ax3 = fig.add_subplot(2,2,4,projection = crs1)

    # Read in the matching Worldview image
    # ------------------------------------
    base_path = '/home/bsorenson/Research/OMI/'
    img_name1 = dt_date_str1.strftime(base_path + 'snapshot-%Y-%m-%dT00_00_00Z_zoomed2.jpg') 
    print("Looking for ",img_name1)

    pic = plt.imread(img_name1)
    ax0.imshow(pic)
    ax0.axis('off')

    base_path = '/home/bsorenson/Research/OMI/'
    img_name2 = dt_date_str2.strftime(base_path + 'snapshot-%Y-%m-%dT00_00_00Z_zoomed2.jpg') 
    print("Looking for ",img_name2)

    pic = plt.imread(img_name2)
    ax3.imshow(pic)
    ax3.axis('off')

    # ------------------------------------------------------------------------
    #
    # Step 3: Plot OMI and CERES data in first two panels
    #
    # ------------------------------------------------------------------------
    lat_lims1 = aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][\
                dt_date_str1.strftime('%H%M')]['Lat']
    lon_lims1 = aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][\
                dt_date_str1.strftime('%H%M')]['Lon']
    lat_lims2 = aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][\
                dt_date_str2.strftime('%H%M')]['Lat']
    lon_lims2 = aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][\
                dt_date_str2.strftime('%H%M')]['Lon']
    ##!#if(zoom):
    ##!#    OMI_base['UVAI'] = np.ma.masked_where(~((OMI_base['LAT'] >= lat_lims[0]) & \
    ##!#                                            (OMI_base['LAT'] <  lat_lims[1]) & \
    ##!#                                            (OMI_base['LON'] >= lon_lims[0]) & \
    ##!#                                            (OMI_base['LON'] <  lon_lims[1])), OMI_base['UVAI'])
    ##!#    #OMI_base['UVAI'] = np.ma.masked_where(OMI_base['UVAI'] > 1, OMI_base['UVAI'])
    ##!#    CERES_hrly[ceres_dtype] = np.ma.masked_where(~((CERES_hrly['lat'] >= lat_lims[0]) & \
    ##!#                                              (CERES_hrly['lat'] <  lat_lims[1]) & \
    ##!#                                              (CERES_hrly['lon'] >= lon_lims[0]) & \
    ##!#                                              (CERES_hrly['lon'] <  lon_lims[1])), CERES_hrly[ceres_dtype])

    # Use the single-swath plotting function to plot OMI data
    # -------------------------------------------------------
    plotOMI_single_swath(ax1, OMI_base1, title = None, \
        circle_bound = circle_bound, gridlines = False)
    plotOMI_single_swath(ax4, OMI_base2, title = None, \
        circle_bound = circle_bound, gridlines = False)

    # Use the single-swath plotting function to plot CERES data
    # ---------------------------------------------------------
    plotCERES_hrly(ax2, CERES_hrly1, 'swf', minlat = minlat, \
        vmin = 90, vmax = 250, title = '', label = 'TOA Flux [Wm$^{-2}$]', \
        circle_bound = circle_bound, gridlines = False, grid_data = True)
    plotCERES_hrly(ax5, CERES_hrly2, 'swf', minlat = minlat, \
        vmin = 90, vmax = 270, title = '', label = 'TOA Flux [Wm$^{-2}$]', \
        circle_bound = circle_bound, gridlines = False, grid_data = True)

    if(zoom):
        lat_lims1 = aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][\
                    dt_date_str1.strftime('%H%M')]['Lat']
        lon_lims1 = aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][\
                    dt_date_str1.strftime('%H%M')]['Lon']
        lat_lims2 = aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][\
                    dt_date_str2.strftime('%H%M')]['Lat']
        lon_lims2 = aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][\
                    dt_date_str2.strftime('%H%M')]['Lon']

        ax1.set_extent([lon_lims1[0], lon_lims1[1], lat_lims1[0], lat_lims1[1]], datacrs)
        ax2.set_extent([lon_lims1[0], lon_lims1[1], lat_lims1[0], lat_lims1[1]], datacrs)
        ax4.set_extent([lon_lims2[0], lon_lims2[1], lat_lims2[0], lat_lims2[1]], datacrs)
        ax5.set_extent([lon_lims2[0], lon_lims2[1], lat_lims2[0], lat_lims2[1]], datacrs)
    else: 
        ax1.set_extent([-180,180,minlat,90], datacrs)
        ax2.set_extent([-180,180,minlat,90], datacrs)
        ax4.set_extent([-180,180,minlat,90], datacrs)
        ax5.set_extent([-180,180,minlat,90], datacrs)

    # If the user wants circles along latitude lines,    
    # plot them here      
    # ----------------------------------------------------
    if(lat_circles is not None):
        plot_lat_circles(ax0, lat_circles) 
        plot_lat_circles(ax1, lat_circles) 

    fig.tight_layout()

    labelsize = 13
    plot_figure_text(ax0, 'MODIS True Color', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax1, 'OMI AI', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'CERES SWF', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'MODIS True Color', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'OMI AI', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'CERES SWF', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = labelsize, backgroundcolor = 'white', halign = 'right')

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax0, '(a)', fontsize = labelsize, backgroundcolor = 'white')
    plot_subplot_label(ax1, '(b)', fontsize = labelsize, backgroundcolor = 'white')
    plot_subplot_label(ax2, '(c)', fontsize = labelsize, backgroundcolor = 'white')
    plot_subplot_label(ax3, '(d)', fontsize = labelsize, backgroundcolor = 'white')
    plot_subplot_label(ax4, '(e)', fontsize = labelsize, backgroundcolor = 'white')
    plot_subplot_label(ax5, '(f)', fontsize = labelsize, backgroundcolor = 'white')

    fig.text(0.015, 0.75, '--- 04:45 UTC 11 Aug 2019 ---', ha='center', va='center', \
        rotation='vertical', weight = 'bold', fontsize = 12)
    fig.text(0.015, 0.26, '--- 18:10 UTC 19 Aug 2017 ---', ha='center', va='center', \
        rotation='vertical', weight = 'bold', fontsize = 12)

    if(save == True):
        outname = 'omi_ceres_compare_'+OMI_date1+'_'+OMI_date2+'_6panel.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

    #return OMI_base, CERES_hrly

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# CERES/NSIDC/NAAPS gridded comparison work
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def plot_compare_CERES_NSIDC_NAAPS_OMI(CERES_data, NSIDC_data, NAAPS_data, \
        OMI_data, \
        month_idx = 4, minlat = 70., trend_type = 'standard', \
        save = False):

    month_obj = datetime(year = 1,month = int(CERES_data['dates'][month_idx][4:]),day = 1)
    str_month = month_obj.strftime('%B')

    # Set up the figure
    #plt.close('all')
    fig1 = plt.figure(figsize = (16,12))
    ax1  = fig1.add_subplot(3,4,5, projection = mapcrs)
    ax2  = fig1.add_subplot(3,4,2, projection = mapcrs)
    ax3  = fig1.add_subplot(3,4,3, projection = mapcrs)
    ax4  = fig1.add_subplot(3,4,4, projection = mapcrs)
    ax5  = fig1.add_subplot(3,4,6)
    ax6  = fig1.add_subplot(3,4,7)
    ax7  = fig1.add_subplot(3,4,8)
    ax8  = fig1.add_subplot(3,4,9, projection = mapcrs)
    ax9  = fig1.add_subplot(3,4,10, projection = mapcrs)
    ax10 = fig1.add_subplot(3,4,11, projection = mapcrs)
    ax11 = fig1.add_subplot(3,4,12, projection = mapcrs)

    cbar_switch = True
    colorbar_label_size = 10
    ceres_trends = plotCERES_MonthTrend(CERES_data,month_idx=month_idx,save=False,\
        trend_type=trend_type,season='sunlight',minlat=minlat,return_trend=True, \
        pax = ax1)
    nsidc_trends = plotNSIDC_MonthTrend(NSIDC_data,month_idx=month_idx,save=False,\
        trend_type=trend_type,season='sunlight',minlat=minlat,\
        return_trend=True, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax2, show_pval = False, uncert_ax = None, title = '')
    naaps_trends = plotNAAPS_MonthTrend(NAAPS_data,month_idx=month_idx,trend_type=trend_type,label = ' ',\
        minlat=minlat,title = ' ', pax = ax3, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, show_pval = False, \
        uncert_ax = None, vmax = 1, return_trend = True)
    omi_trends = plotOMI_MonthTrend(OMI_data,month_idx=month_idx,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        return_trend = True, \
        minlat=minlat,title = None, pax = ax4)

    if(np.max(CERES_data['lon']) > 270):

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
            ceres_trends[ii,:] = np.concatenate([ceres_trends[ii,:][over_180],\
                ceres_trends[ii,:][under_180]])
   

    mask_nsidc = np.ma.masked_invalid(nsidc_trends)
    mask_nsidc = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_nsidc)
    mask_naaps = np.ma.masked_invalid(naaps_trends)
    mask_naaps = np.ma.masked_where(mask_nsidc == 0., mask_naaps)
    mask_omi   = np.ma.masked_invalid(omi_trends)
    mask_omi   = np.ma.masked_where(mask_nsidc == 0., mask_omi)
    mask_ceres = np.ma.masked_invalid(ceres_trends)
    mask_ceres = np.ma.masked_where(mask_nsidc == 0., mask_ceres)

    ax8.pcolormesh(local_lon, local_lat, mask_ceres, \
        transform = datacrs, vmax = 30, vmin = -30, shading = 'auto', \
        cmap = 'bwr')
    ax8.coastlines()
    ax8.set_extent([-180,180,minlat, 90], datacrs)
    ax9.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        mask_nsidc, cmap = 'bwr', \
        transform = datacrs, vmax = 30, vmin = -30, shading = 'auto')
    ax9.coastlines()
    ax9.set_extent([-180,180,minlat, 90], datacrs)
    ax10.pcolormesh(NAAPS_data['lons'], NAAPS_data['lats'], \
        mask_naaps, cmap = 'bwr', \
        transform = datacrs, vmax = 1, vmin = -1, shading = 'auto')
    ax10.coastlines()
    ax10.set_extent([-180,180,minlat, 90], datacrs)
    ax11.pcolormesh(OMI_data['LON'], OMI_data['LAT'], \
        mask_omi, cmap = 'bwr', \
        transform = datacrs, vmax = 0.5, vmin = -0.5, shading = 'auto')
    ax11.coastlines()
    ax11.set_extent([-180,180,minlat, 90], datacrs)

    final_ceres = ceres_trends[(~mask_ceres.mask) & (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & (~mask_omi.mask)]
    final_nsidc = nsidc_trends[(~mask_ceres.mask) & (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & (~mask_omi.mask)]
    final_naaps = naaps_trends[(~mask_ceres.mask) & (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & (~mask_omi.mask)]
    final_omi   = omi_trends[(~mask_ceres.mask) & (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & (~mask_omi.mask)]

    #return ceres_trends,nsidc_trends,naaps_trends

    ax5.scatter(final_nsidc.flatten(), final_ceres.flatten(), s = 6, color = 'k')
    ax6.scatter(final_naaps.flatten(), final_ceres.flatten(), s = 6, color = 'k')
    ax7.scatter(final_omi.flatten(),   final_ceres.flatten(), s = 6, color = 'k')
    ax6.set_xlim(ax6.get_xlim()[0], 0.5)

    ax2.set_title('NSIDC Sea Ice Concentration Trends')
    ax3.set_title('NAAPS Near-Surface\nSmoke Aerosol Trends')
    ax4.set_title('OMI UVAI Trends')
    
    fig1.tight_layout()
    fig1.tight_layout()
    if(save):
        outname = '_'.join(['ceres','nsidc','naaps','omi','trend','comp',\
            month_obj.strftime('%b'), CERES_data['param'].split('_')[2]]) + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_compare_CERES_NSIDC_NAAPS_OMI_combined(CERES_swclr, CERES_swall, \
        NSIDC_data, NAAPS_data, OMI_data, \
        month_idx = 4, minlat = 70., trend_type = 'standard', \
        save = False):

    month_obj = datetime(year = 1,month = int(CERES_data['dates'][month_idx][4:]),day = 1)
    str_month = month_obj.strftime('%B')

    # Set up the figure
    plt.close('all')
    fig1 = plt.figure(figsize = (12,6))
    ax1  = fig1.add_subplot(2,4,1, projection = mapcrs)
    ax2  = fig1.add_subplot(2,4,2, projection = mapcrs)
    ax3  = fig1.add_subplot(2,4,3, projection = mapcrs)
    ax4  = fig1.add_subplot(2,4,4, projection = mapcrs)
    ax5  = fig1.add_subplot(2,4,5, projection = mapcrs)
    ax6  = fig1.add_subplot(2,4,6)
    ax7  = fig1.add_subplot(2,4,7)
    ax8  = fig1.add_subplot(2,4,8)

    ceres_clr_trends, ceres_clr_pvals, ceres_clr_uncert = \
        calcCERES_grid_trend(CERES_swclr, month_idx, trend_type, \
        minlat)
    ceres_all_trends, ceres_all_pvals, ceres_all_uncert = \
        calcCERES_grid_trend(CERES_swall, month_idx, trend_type, \
        minlat)
    naaps_trends, naaps_pvals, naaps_uncert = \
        calcNAAPS_grid_trend(NAAPS_data, month_idx, trend_type, \
        minlat)
    nsidc_trends, nsidc_pvals, nsidc_uncert = \
        calcNSIDC_grid_trend(NSIDC_data, month_idx, trend_type, \
        minlat)
    omi_trends, omi_pvals, omi_uncert = \
        calcOMI_grid_trend(OMI_data, month_idx, trend_type, \
        minlat)
    
    ##!#cbar_switch = True
    ##!#colorbar_label_size = 10
    ##!#ceres_clr_trends = plotCERES_MonthTrend(CERES_swclr,month_idx=month_idx,save=False,\
    ##!#    trend_type=trend_type,season='sunlight',minlat=minlat,return_trend=True, \
    ##!#    pax = ax1)
    ##!#ceres_all_trends = plotCERES_MonthTrend(CERES_swall,month_idx=month_idx,save=False,\
    ##!#    trend_type=trend_type,season='sunlight',minlat=minlat,return_trend=True, \
    ##!#    pax = ax2)
    ##!#nsidc_trends = plotNSIDC_MonthTrend(NSIDC_data,month_idx=month_idx,save=False,\
    ##!#    trend_type=trend_type,season='sunlight',minlat=minlat,\
    ##!#    return_trend=True, colorbar = cbar_switch, \
    ##!#    colorbar_label_size = colorbar_label_size, \
    ##!#    pax = ax3, show_pval = False, uncert_ax = None, title = '')
    ##!#naaps_trends = plotNAAPS_MonthTrend(NAAPS_data,month_idx=month_idx,trend_type=trend_type,label = ' ',\
    ##!#    minlat=minlat,title = ' ', pax = ax4, colorbar = cbar_switch, \
    ##!#    colorbar_label_size = colorbar_label_size, show_pval = False, \
    ##!#    uncert_ax = None, vmax = 1, return_trend = True)
    ##!#omi_trends = plotOMI_MonthTrend(OMI_data,month_idx=month_idx,\
    ##!#    trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
    ##!#    return_trend = True, \
    ##!#    minlat=minlat,title = None, pax = ax5)

    if(np.max(CERES_swclr['lon']) > 270):

        # Flip the CERES data to convert the longitudes from 0 - 360 to -180 - 180
        # ------------------------------------------------------------------------
        local_lon = np.copy(CERES_swclr['lon'])
        local_lat = np.copy(CERES_swclr['lat'])
        over_180  = np.where(CERES_swclr['lon'][0,:] >= 179.9999)
        under_180 = np.where(CERES_swclr['lon'][0,:] < 179.9999)

        for ii in range(local_lon.shape[0]):
            local_lon[ii,:] = np.concatenate([local_lon[ii,:][over_180] - 360.,\
                local_lon[ii,:][under_180]])
            local_lat[ii,:] = np.concatenate([local_lat[ii,:][over_180],\
                local_lat[ii,:][under_180]])
            ceres_clr_trends[ii,:] = np.concatenate([ceres_clr_trends[ii,:][over_180],\
                ceres_clr_trends[ii,:][under_180]])
            ceres_all_trends[ii,:] = np.concatenate([ceres_all_trends[ii,:][over_180],\
                ceres_all_trends[ii,:][under_180]])

    mask_nsidc = np.ma.masked_invalid(nsidc_trends)
    mask_nsidc = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_nsidc)
    mask_naaps = np.ma.masked_invalid(naaps_trends)
    mask_naaps = np.ma.masked_where(naaps_trends > 1, mask_naaps)
    mask_naaps = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_naaps)
    mask_omi   = np.ma.masked_invalid(omi_trends)
    mask_omi   = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_omi)
    mask_ceres_clr = np.ma.masked_invalid(ceres_clr_trends)
    mask_ceres_clr = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_ceres_clr)
    mask_ceres_all = np.ma.masked_invalid(ceres_all_trends)
    mask_ceres_all = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_ceres_all)

    labelsize = 9
    # Map 1: CERES SWCLR
    mesh = ax1.pcolormesh(local_lon, local_lat, mask_ceres_clr, \
        transform = datacrs, vmax = 30, vmin = -30, shading = 'auto', \
        cmap = 'bwr')
    cbar = plt.colorbar(mesh,ax = ax1, orientation='vertical',shrink = 0.8, \
        extend = 'both')
    cbar.set_label('SW Flux Trend [W m$^{-2}$]', \
        fontsize = labelsize)
    ax1.coastlines()
    ax1.set_extent([-180,180,minlat, 90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    # Map 2: CERES SWALL
    mesh = ax2.pcolormesh(local_lon, local_lat, mask_ceres_all, \
        transform = datacrs, vmax = 30, vmin = -30, shading = 'auto', \
        cmap = 'bwr')
    cbar = plt.colorbar(mesh,ax = ax2, orientation='vertical',shrink = 0.8, \
        extend = 'both')
    cbar.set_label('SW Flux Trend [W m$^{-2}$]', \
        fontsize = labelsize)
    ax2.coastlines()
    ax2.set_extent([-180,180,minlat, 90], datacrs)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    # Map 3: NSIDC
    mesh = ax3.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        mask_nsidc, cmap = 'bwr', \
        transform = datacrs, vmax = 30, vmin = -30, shading = 'auto')
    cbar = plt.colorbar(mesh,ax = ax3, orientation='vertical',shrink = 0.8, \
        extend = 'both')
    cbar.set_label('Sea Ice Trend [%]', \
        fontsize = labelsize)
    ax3.coastlines()
    ax3.set_extent([-180,180,minlat, 90], datacrs)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    # Map 4: NAAPS
    mesh = ax4.pcolormesh(NAAPS_data['lons'], NAAPS_data['lats'], \
        mask_naaps, cmap = 'bwr', \
        transform = datacrs, vmax = 0.5, vmin = -0.5, shading = 'auto')
    cbar = plt.colorbar(mesh,ax = ax4, orientation='vertical',shrink = 0.8, \
        extend = 'both')
    cbar.set_label('Smoke Trend [μg m$^{-3}$]', \
        fontsize = labelsize)
    ax4.coastlines()
    ax4.set_extent([-180,180,minlat, 90], datacrs)
    ax4.set_boundary(circle, transform=ax4.transAxes)
    # Map 5: OMI
    mesh = ax5.pcolormesh(OMI_data['LON'], OMI_data['LAT'], \
        mask_omi, cmap = 'bwr', \
        transform = datacrs, vmax = 0.5, vmin = -0.5, shading = 'auto')
    cbar = plt.colorbar(mesh,ax = ax5, orientation='vertical',shrink = 0.8, \
        extend = 'both')
    cbar.set_label('UVAI Trend [AI]', \
        fontsize = labelsize)
    ax5.coastlines()
    ax5.set_extent([-180,180,minlat, 90], datacrs)
    ax5.set_boundary(circle, transform=ax5.transAxes)

    titlesize = 9
    ax1.set_title('Aqua CERES\nClear-sky SWF Trends', fontsize = titlesize)
    ax2.set_title('Aqua CERES\nAll-sky SWF Trends', fontsize = titlesize)
    ax3.set_title('SSMI/S Sea Ice\nConcentration Trends', fontsize = titlesize)
    ax4.set_title('NAAPS-RA Near-surface\nSmoke Aerosol Trends', fontsize = titlesize)
    ax5.set_title('OMI Perturbed\nUVAI Trends', fontsize = titlesize)

    add_gridlines(ax1)
    add_gridlines(ax2)
    add_gridlines(ax3)
    add_gridlines(ax4)
    add_gridlines(ax5)

    final_ceres_clr = ceres_clr_trends[(~mask_ceres_clr.mask) & \
        (~mask_ceres_all.mask) & \
        (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & \
        (~mask_omi.mask)]
    final_ceres_all = ceres_all_trends[(~mask_ceres_clr.mask) & \
        (~mask_ceres_all.mask) & \
        (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & \
        (~mask_omi.mask)]
    final_nsidc = nsidc_trends[(~mask_ceres_clr.mask) & \
        (~mask_ceres_all.mask) & \
        (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & \
        (~mask_omi.mask)]
    final_naaps = naaps_trends[(~mask_ceres_clr.mask) & \
        (~mask_ceres_all.mask) & \
        (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & \
        (~mask_omi.mask)]
    final_omi   = omi_trends[(~mask_ceres_clr.mask) & \
        (~mask_ceres_all.mask) & \
        (~mask_nsidc.mask) & \
        (~mask_naaps.mask) & \
        (~mask_omi.mask)]

    #return ceres_trends,nsidc_trends,naaps_trends

    # Scatter 1: NSIDC vs CERES SWCLR
    sizer = 4
    ax6.scatter(final_nsidc.flatten(), final_ceres_clr.flatten(), s = sizer, color = 'k')
    # Scatter 2: OMI vs CERES SWALL
    ax7.scatter(final_omi.flatten(), final_ceres_all.flatten(), s = sizer, color = 'k')
    # Scatter 3: NSIDC vs CERES SWCLR
    ax8.scatter(final_naaps.flatten(),  final_ceres_clr.flatten(), s = sizer, color = 'k')
    #ax8.set_xlim(ax8.get_xlim()[0], 0.5)

    ax6_result = stats.linregress(final_nsidc.flatten(), final_ceres_clr.flatten())
    ax7_result = stats.linregress(final_omi.flatten(), final_ceres_all.flatten())
    ax8_result = stats.linregress(final_naaps.flatten(), final_ceres_clr.flatten())
    
    print(ax6_result)

    corr1 = ax6_result.rvalue**2.
    corr2 = ax7_result.rvalue**2.
    corr3 = ax8_result.rvalue**2.

    #ax6.text(0,    80, 'r$^{2}$ = ' + str(np.round(corr1, 3)))
    #ax7.text(0.05, 35, 'r$^{2}$ = ' + str(np.round(corr2, 3)), backgroundcolor = 'white')
    #ax8.text(0.5,  80, 'r$^{2}$ = ' + str(np.round(corr3, 3)))

    plot_figure_text(ax6, 'r$^{2}$ = ' + str(np.round(corr1, 3)), xval = None, yval = None, transform = None, \
        color = 'black', fontsize = 10, backgroundcolor = 'white',\
        halign = 'right', location = 'lower_right')
    plot_figure_text(ax7, 'r$^{2}$ = ' + str(np.round(corr2, 3)), xval = None, yval = None, transform = None, \
        color = 'black', fontsize = 10, backgroundcolor = 'white',\
        halign = 'right', location = 'upper_right')
    plot_figure_text(ax8, 'r$^{2}$ = ' + str(np.round(corr3, 3)), xval = None, yval = None, transform = None, \
        color = 'black', fontsize = 10, backgroundcolor = 'white',\
        halign = 'right', location = 'upper_right')

    print("Correlation 1", ax6_result.rvalue**2.)
    print("Correlation 2", ax7_result.rvalue**2.)
    print("Correlation 3", ax8_result.rvalue**2.)

    sloper = 'linregress'
    plot_trend_line(ax6, final_nsidc.flatten(), final_ceres_clr.flatten(), \
        color='r', linestyle = '-', slope = sloper)
    plot_trend_line(ax7, final_omi.flatten(), final_ceres_all.flatten(), \
        color='r', linestyle = '-', slope = sloper)
    plot_trend_line(ax8, final_naaps.flatten(), final_ceres_clr.flatten(), \
        color='r', linestyle = '-', slope = sloper)

    ax6.set_title('Sea Ice Concentration Trend vs.\nCERES Clear-sky SWF Trend', fontsize = titlesize)
    ax7.set_title('OMI Perturbed UVAI Trend vs.\nCERES All-sky SWF Trend', fontsize = titlesize)
    ax8.set_title('NAAPS-RA Surface Smoke Trend vs.\nCERES Clear-sky SWF Trend', fontsize = titlesize)

    ax6.set_xlabel('Sea Ice Trend [% / Study Period]', fontsize = labelsize)
    ax7.set_xlabel('UVAI Trend [AI/Study Period]', fontsize = labelsize)
    ax8.set_xlabel('Near-surface Smoke Trend [μg m$^{-3}$]', fontsize = labelsize)
    ax6.set_ylabel('Clear-sky SWF Trend [W m$^{-2}$]', fontsize = labelsize)
    ax7.set_ylabel('All-sky SWF Trend [W m$^{-2}$]', fontsize = labelsize)
    ax8.set_ylabel('Clear-sky SWF Trend [W m$^{-2}$]', fontsize = labelsize)

    plot_subplot_label(ax1, '(a)', location = 'upper_left')
    plot_subplot_label(ax2, '(b)', location = 'upper_left')
    plot_subplot_label(ax3, '(c)', location = 'upper_left')
    plot_subplot_label(ax4, '(d)', location = 'upper_left')
    plot_subplot_label(ax5, '(e)', location = 'upper_left')
    plot_subplot_label(ax6, '(f)', location = 'upper_left')
    plot_subplot_label(ax7, '(g)', location = 'upper_left')
    plot_subplot_label(ax8, '(h)', location = 'upper_left')
   
    plt.suptitle(month_obj.strftime("%B Monthly Trends (2005 - 2020)"))
 
    fig1.tight_layout()
    if(save):
        outname = '_'.join(['ceres','nsidc','naaps','omi','trend','comp',\
            month_obj.strftime('%b'),'combined']) + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# Assumed that the NAAPS_data dictionary is returned from 
# readgridNAAPS_monthly
def plot_compare_CERES_NAAPS_monthly(CERES_swclr, CERES_swall, \
        NSIDC_data, NAAPS_data, param = 'smoke_conc_sfc', \
        month_idx = 4, minlat = 70., trend_type = 'standard', \
        save = False):

    month_obj = datetime(year = 1,month = int(CERES_data['dates'][month_idx][4:]),day = 1)
    str_month = month_obj.strftime('%B')

    # Set up the figure
    plt.close('all')
    fig1 = plt.figure(figsize = (9,6))
    ax6  = fig1.add_subplot(2,3,1, projection = mapcrs) # CERES All-sky
    ax1  = fig1.add_subplot(2,3,3, projection = mapcrs) # CERES All-sky
    ax2  = fig1.add_subplot(2,3,2, projection = mapcrs) # CERES clear-sky
    ax3  = fig1.add_subplot(2,3,4, projection = mapcrs) # NAAPS variable
    ax4  = fig1.add_subplot(2,3,6)
    ax5  = fig1.add_subplot(2,3,5)
    #fig1 = plt.figure(figsize = (12,6))
    #ax1  = fig1.add_subplot(2,4,1, projection = mapcrs)
    #ax2  = fig1.add_subplot(2,4,2, projection = mapcrs)
    #ax3  = fig1.add_subplot(2,4,3, projection = mapcrs)
    #ax4  = fig1.add_subplot(2,4,4, projection = mapcrs)
    #ax5  = fig1.add_subplot(2,4,5, projection = mapcrs)
    #ax6  = fig1.add_subplot(2,4,6)
    #ax7  = fig1.add_subplot(2,4,7)
    #ax8  = fig1.add_subplot(2,4,8)

    ceres_clr_trends, ceres_clr_pvals, ceres_clr_uncert = \
        calcCERES_grid_trend(CERES_swclr, month_idx, trend_type, \
        minlat)
    ceres_all_trends, ceres_all_pvals, ceres_all_uncert = \
        calcCERES_grid_trend(CERES_swall, month_idx, trend_type, \
        minlat)
    nsidc_trends, nsidc_pvals, nsidc_uncert = \
        calcNSIDC_grid_trend(NSIDC_data, month_idx, trend_type, \
        minlat)
    naaps_trends, naaps_pvals, naaps_uncert = \
        calcNAAPS_grid_trend(NAAPS_data, month_idx, trend_type, \
        minlat, param = param)


    #nsidc_trends, nsidc_pvals, nsidc_uncert = \
    #    calcNSIDC_grid_trend(NSIDC_data, month_idx, trend_type, \
    #    minlat)
    #omi_trends, omi_pvals, omi_uncert = \
    #    calcOMI_grid_trend(OMI_data, month_idx, trend_type, \
    #    minlat)
    
    ##!#cbar_switch = True
    ##!#colorbar_label_size = 10
    ##!#ceres_clr_trends = plotCERES_MonthTrend(CERES_swclr,month_idx=month_idx,save=False,\
    ##!#    trend_type=trend_type,season='sunlight',minlat=minlat,return_trend=True, \
    ##!#    pax = ax1)
    ##!#ceres_all_trends = plotCERES_MonthTrend(CERES_swall,month_idx=month_idx,save=False,\
    ##!#    trend_type=trend_type,season='sunlight',minlat=minlat,return_trend=True, \
    ##!#    pax = ax2)
    ##!#nsidc_trends = plotNSIDC_MonthTrend(NSIDC_data,month_idx=month_idx,save=False,\
    ##!#    trend_type=trend_type,season='sunlight',minlat=minlat,\
    ##!#    return_trend=True, colorbar = cbar_switch, \
    ##!#    colorbar_label_size = colorbar_label_size, \
    ##!#    pax = ax3, show_pval = False, uncert_ax = None, title = '')
    ##!#naaps_trends = plotNAAPS_MonthTrend(NAAPS_data,month_idx=month_idx,trend_type=trend_type,label = ' ',\
    ##!#    minlat=minlat,title = ' ', pax = ax4, colorbar = cbar_switch, \
    ##!#    colorbar_label_size = colorbar_label_size, show_pval = False, \
    ##!#    uncert_ax = None, vmax = 1, return_trend = True)
    ##!#omi_trends = plotOMI_MonthTrend(OMI_data,month_idx=month_idx,\
    ##!#    trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
    ##!#    return_trend = True, \
    ##!#    minlat=minlat,title = None, pax = ax5)

    if(np.max(CERES_swclr['lon']) > 270):

        # Flip the CERES data to convert the longitudes from 0 - 360 to -180 - 180
        # ------------------------------------------------------------------------
        local_lon = np.copy(CERES_swclr['lon'])
        local_lat = np.copy(CERES_swclr['lat'])
        over_180  = np.where(CERES_swclr['lon'][0,:] >= 179.9999)
        under_180 = np.where(CERES_swclr['lon'][0,:] < 179.9999)

        for ii in range(local_lon.shape[0]):
            local_lon[ii,:] = np.concatenate([local_lon[ii,:][over_180] - 360.,\
                local_lon[ii,:][under_180]])
            local_lat[ii,:] = np.concatenate([local_lat[ii,:][over_180],\
                local_lat[ii,:][under_180]])
            ceres_clr_trends[ii,:] = np.concatenate([ceres_clr_trends[ii,:][over_180],\
                ceres_clr_trends[ii,:][under_180]])
            ceres_all_trends[ii,:] = np.concatenate([ceres_all_trends[ii,:][over_180],\
                ceres_all_trends[ii,:][under_180]])

    mask_nsidc = np.ma.masked_invalid(nsidc_trends)
    mask_nsidc = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_nsidc)
    mask_naaps = np.ma.masked_invalid(naaps_trends)
    #mask_naaps = np.ma.masked_where(naaps_trends > 1, mask_naaps)
    mask_naaps = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_naaps)
    #mask_omi   = np.ma.masked_invalid(omi_trends)
    #mask_omi   = np.ma.masked_where((mask_nsidc == 0.) | \
    #    (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_omi)
    mask_ceres_clr = np.ma.masked_invalid(ceres_clr_trends)
    mask_ceres_clr = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_ceres_clr)
    mask_ceres_all = np.ma.masked_invalid(ceres_all_trends)
    mask_ceres_all = np.ma.masked_where((mask_nsidc == 0.) | \
        (NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40), mask_ceres_all)

    #plotNSIDC_MonthClimo(NSIDC_data,month_idx,minlat = minlat, ax = ax6, title = '', \
    #    colorbar = True)

    mask_data = np.ma.masked_where(\
        NSIDC_data['MONTH_CLIMO'][month_idx,:,:] < 40, 
        NSIDC_data['MONTH_CLIMO'][month_idx,:,:])
    labelsize = 9
    # Map 1: NSIDC ICE
    mesh = ax6.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
        mask_data, \
        transform = datacrs, vmax = 10, vmin = 100, shading = 'auto', \
        cmap = 'ocean')
    cbar = plt.colorbar(mesh,ax = ax6, orientation='vertical', \
        shrink = 0.8, \
        extend = 'both')
    cbar.set_label('Sea Ice Concentration [%]', \
        fontsize = labelsize)
    ax6.coastlines()
    ax6.set_extent([-180,180,minlat, 90], datacrs)
    ax6.set_boundary(circle, transform=ax6.transAxes)

    # Map 1: CERES SWCLR
    mesh = ax1.pcolormesh(local_lon, local_lat, mask_ceres_clr, \
        transform = datacrs, vmax = 30, vmin = -30, shading = 'auto', \
        cmap = 'bwr')
    cbar = plt.colorbar(mesh,ax = ax1, orientation='vertical', \
        shrink = 0.8, \
        extend = 'both')
    cbar.set_label('SW Flux Trend [W m$^{-2}$]', \
        fontsize = labelsize)
    ax1.coastlines()
    ax1.set_extent([-180,180,minlat, 90], datacrs)
    ax1.set_boundary(circle, transform=ax1.transAxes)
    # Map 2: CERES SWALL
    mesh = ax2.pcolormesh(local_lon, local_lat, mask_ceres_all, \
        transform = datacrs, vmax = 30, vmin = -30, shading = 'auto', \
        cmap = 'bwr')
    cbar = plt.colorbar(mesh,ax = ax2, orientation='vertical', \
        shrink = 0.8, \
        extend = 'both')
    cbar.set_label('SW Flux Trend [W m$^{-2}$]', \
        fontsize = labelsize)
    ax2.coastlines()
    ax2.set_extent([-180,180,minlat, 90], datacrs)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ## Map 3: NSIDC
    #mesh = ax3.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    #    mask_nsidc, cmap = 'bwr', \
    #    transform = datacrs, vmax = 30, vmin = -30, shading = 'auto')
    #cbar = plt.colorbar(mesh,ax = ax3, orientation='vertical',shrink = 0.8, \
    #    extend = 'both')
    #cbar.set_label('Sea Ice Trend [%]', \
    #    fontsize = labelsize)
    #ax3.coastlines()
    #ax3.set_extent([-180,180,minlat, 90], datacrs)
    #ax3.set_boundary(circle, transform=ax3.transAxes)
    # Map 4: NAAPS

    lim_dict = {
        'smoke_conc_sfc': 1.5,
        'smoke_conc_sfc_yearly': 5., 
        'smoke_totsink': 1.5,
        'smoke_totsink_yearly': 3.,
        'smoke_drysink': 0.1,
        'smoke_drysink_yearly': 0.25,
        'smoke_wetsink': 1,
        'smoke_wetsink_yearly': 3.,
    }   

    vmin = -lim_dict[param]
    vmax = lim_dict[param]

    mask_naaps = np.ma.masked_where((mask_naaps > 3*vmax) | \
        (mask_naaps < 3*vmin), mask_naaps)

    mesh = ax3.pcolormesh(NAAPS_data['lons'], NAAPS_data['lats'], \
        mask_naaps, cmap = 'bwr', \
        transform = datacrs, vmax = vmax, vmin = vmin, shading = 'auto')
    cbar = plt.colorbar(mesh,ax = ax3, orientation='vertical', \
        shrink = 0.8, \
        extend = 'both')
    cbar.set_label('Smoke Trend [μg m$^{-3}$]', \
        fontsize = labelsize)
    ax3.coastlines()
    ax3.set_extent([-180,180,minlat, 90], datacrs)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ## Map 5: OMI
    #mesh = ax5.pcolormesh(OMI_data['LON'], OMI_data['LAT'], \
    #    mask_omi, cmap = 'bwr', \
    #    transform = datacrs, vmax = 0.5, vmin = -0.5, shading = 'auto')
    #cbar = plt.colorbar(mesh,ax = ax5, orientation='vertical',shrink = 0.8, \
    #    extend = 'both')
    #cbar.set_label('UVAI Trend [AI]', \
    #    fontsize = labelsize)
    #ax5.coastlines()
    #ax5.set_extent([-180,180,minlat, 90], datacrs)
    #ax5.set_boundary(circle, transform=ax5.transAxes)

    titlesize = 9
    ax6.set_title('SSMI/S Average\nSea Ice Conc.', fontsize = titlesize)
    ax1.set_title('Aqua CERES\nClear-sky SWF Trends', fontsize = titlesize)
    ax2.set_title('Aqua CERES\nAll-sky SWF Trends', fontsize = titlesize)
    #ax3.set_title('SSMI/S Sea Ice\nConcentration Trends', fontsize = titlesize)
    ax3.set_title('NAAPS-RA Near-surface\n' + param, fontsize = titlesize)
    #ax5.set_title('OMI Perturbed\nUVAI Trends', fontsize = titlesize)

    add_gridlines(ax1)
    add_gridlines(ax2)
    add_gridlines(ax3)
    #add_gridlines(ax4)
    #add_gridlines(ax5)

    final_ceres_clr = ceres_clr_trends[(~mask_ceres_clr.mask) & \
        (~mask_ceres_all.mask) & \
        (~mask_naaps.mask)]
    final_ceres_all = ceres_all_trends[(~mask_ceres_clr.mask) & \
        (~mask_ceres_all.mask) & \
        (~mask_naaps.mask)]
    final_naaps = naaps_trends[(~mask_ceres_clr.mask) & \
        (~mask_ceres_all.mask) & \
        (~mask_naaps.mask)]
    #final_ceres_clr = ceres_clr_trends[(~mask_ceres_clr.mask) & \
    #    (~mask_ceres_all.mask) & \
    #    (~mask_nsidc.mask) & \
    #    (~mask_naaps.mask) & \
    #    (~mask_omi.mask)]
    #final_ceres_all = ceres_all_trends[(~mask_ceres_clr.mask) & \
    #    (~mask_ceres_all.mask) & \
    #    (~mask_nsidc.mask) & \
    #    (~mask_naaps.mask) & \
    #    (~mask_omi.mask)]
    #final_nsidc = nsidc_trends[(~mask_ceres_clr.mask) & \
    #    (~mask_ceres_all.mask) & \
    #    (~mask_nsidc.mask) & \
    #    (~mask_naaps.mask) & \
    #    (~mask_omi.mask)]
    #final_naaps = naaps_trends[(~mask_ceres_clr.mask) & \
    #    (~mask_ceres_all.mask) & \
    #    (~mask_nsidc.mask) & \
    #    (~mask_naaps.mask) & \
    #    (~mask_omi.mask)]
    #final_omi   = omi_trends[(~mask_ceres_clr.mask) & \
    #    (~mask_ceres_all.mask) & \
    #    (~mask_nsidc.mask) & \
    #    (~mask_naaps.mask) & \
    #    (~mask_omi.mask)]

    #return ceres_trends,nsidc_trends,naaps_trends

    # Scatter 1: NSIDC vs CERES SWCLR
    sizer = 4
    ax4.scatter(final_naaps.flatten(), final_ceres_clr.flatten(), s = sizer, color = 'k')
    # Scatter 2: OMI vs CERES SWALL
    ax5.scatter(final_naaps.flatten(), final_ceres_all.flatten(), s = sizer, color = 'k')
    # Scatter 3: NSIDC vs CERES SWCLR
    #ax8.scatter(final_naaps.flatten(),  final_ceres_clr.flatten(), s = sizer, color = 'k')
    #ax8.set_xlim(ax8.get_xlim()[0], 0.5)

    ax4_result = stats.linregress(final_naaps.flatten(), final_ceres_clr.flatten())
    ax5_result = stats.linregress(final_naaps.flatten(), final_ceres_all.flatten())
    #ax8_result = stats.linregress(final_naaps.flatten(), final_ceres_clr.flatten())
    
    print(ax4_result)
    print(ax5_result)

    corr1 = ax4_result.rvalue**2.
    corr2 = ax5_result.rvalue**2.

    #ax6.text(0,    80, 'r$^{2}$ = ' + str(np.round(corr1, 3)))
    #ax7.text(0.05, 35, 'r$^{2}$ = ' + str(np.round(corr2, 3)), backgroundcolor = 'white')
    #ax8.text(0.5,  80, 'r$^{2}$ = ' + str(np.round(corr3, 3)))

    plot_figure_text(ax4, 'r$^{2}$ = ' + str(np.round(corr1, 3)), xval = None, yval = None, transform = None, \
        color = 'black', fontsize = 10, backgroundcolor = 'white',\
        halign = 'right', location = 'lower_right')
    plot_figure_text(ax5, 'r$^{2}$ = ' + str(np.round(corr2, 3)), xval = None, yval = None, transform = None, \
        color = 'black', fontsize = 10, backgroundcolor = 'white',\
        halign = 'right', location = 'upper_right')
    #plot_figure_text(ax8, 'r$^{2}$ = ' + str(np.round(corr3, 3)), xval = None, yval = None, transform = None, \
    #    color = 'black', fontsize = 10, backgroundcolor = 'white',\
    #    halign = 'right', location = 'upper_right')

    print("Correlation 1", ax4_result.rvalue**2.)
    print("Correlation 2", ax5_result.rvalue**2.)

    sloper = trend_type
    plot_trend_line(ax4, final_naaps.flatten(), final_ceres_clr.flatten(), \
        color='r', linestyle = '-', slope = sloper)
    plot_trend_line(ax5, final_naaps.flatten(), final_ceres_all.flatten(), \
        color='r', linestyle = '-', slope = sloper)
    #plot_trend_line(ax8, final_naaps.flatten(), final_ceres_clr.flatten(), \
    #    color='r', linestyle = '-', slope = sloper)

    ax4.set_title('NAAPS ' + param + ' Trend vs.\nCERES Clear-sky SWF Trend', fontsize = titlesize)
    ax5.set_title('NAAPS ' + param + ' Trend vs.\nCERES All-sky SWF Trend', fontsize = titlesize)
    #ax8.set_title('NAAPS-RA Surface Smoke Trend vs.\nCERES Clear-sky SWF Trend', fontsize = titlesize)

    ax4.set_xlabel(param + ' Trend', fontsize = labelsize)
    ax5.set_xlabel(param + ' Trend', fontsize = labelsize)
    ax4.set_ylabel('Clear-sky SWF Trend [W m$^{-2}$]', fontsize = labelsize)
    ax5.set_ylabel('All-sky SWF Trend [W m$^{-2}$]', fontsize = labelsize)

    plot_subplot_label(ax6, '(a)', location = 'upper_left')
    plot_subplot_label(ax1, '(b)', location = 'upper_left')
    plot_subplot_label(ax2, '(c)', location = 'upper_left')
    plot_subplot_label(ax3, '(d)', location = 'upper_left')
    plot_subplot_label(ax4, '(e)', location = 'upper_left')
    plot_subplot_label(ax5, '(f)', location = 'upper_left')
    #plot_subplot_label(ax6, '(f)', location = 'upper_left')
    #plot_subplot_label(ax7, '(g)', location = 'upper_left')
    #plot_subplot_label(ax8, '(h)', location = 'upper_left')
   
    plt.suptitle(month_obj.strftime("%B Monthly Trends (2005 - 2019)"))
 
    fig1.tight_layout()
    if(save):
        outname = '_'.join(['ceres','naaps',param, 'monthly','trend','comp',\
            month_obj.strftime('%b'),'combined']) + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# Assumed that the NAAPS_data dictionary is returned from 
# readgridNAAPS_monthly
def plot_compare_CERES_NAAPS_monthly_timeseries(CERES_swclr, CERES_swall, \
        NSIDC_data, NAAPS_data, param = 'smoke_conc_sfc', \
        month_idx = 4, minlat = 70.5, maxlat = 80.5, trend_type = 'standard', \
        save = False):

    month_obj = datetime(year = 1,month = int(CERES_data['dates'][month_idx][4:]),day = 1)
    str_month = month_obj.strftime('%B')

    # Prep the data
    naaps_idx = np.where((NAAPS_data['lats'][:,0] >= minlat) & (NAAPS_data['lats'][:,0] <= maxlat))
    ceres_idx = np.where((CERES_swclr['lat'][:,0] >= minlat) & (CERES_swclr['lat'][:,0] <= maxlat))
    nsidc_idx = np.where((NSIDC_data['grid_lat'][:,0] >= minlat) & (NSIDC_data['grid_lat'][:,0] <= maxlat))
  
    local_CERES_swclr = np.ma.masked_where(CERES_swclr['data'] == -999., CERES_swclr['data'])
    local_CERES_swall = np.ma.masked_where(CERES_swall['data'] == -999., CERES_swall['data'])
 
    # Average the data along latitudes
    # -------------------------------- 
    zonal_naaps = np.mean(NAAPS_data[param][4::6,naaps_idx[0],:], axis = 1)
    zonal_ceres_swclr = np.mean(local_CERES_swclr[4::6,naaps_idx[0],:], axis = 1)
    zonal_ceres_swall = np.mean(local_CERES_swall[4::6,naaps_idx[0],:], axis = 1)
    zonal_nsidc = np.mean(NSIDC_data['grid_ice_conc'][4::6,naaps_idx[0],:], axis = 1)

    ## Adjust the CERES longitudes before averaging atain
    #if(np.max(CERES_swclr['lon']) > 270):

    #    # Flip the CERES data to convert the longitudes from 0 - 360 to -180 - 180
    #    # ------------------------------------------------------------------------
    #    local_lon = np.copy(CERES_swclr['lon'])
    #    local_lat = np.copy(CERES_swclr['lat'])
    #    over_180  = np.where(CERES_swclr['lon'][0,:] >= 179.9999)
    #    under_180 = np.where(CERES_swclr['lon'][0,:] < 179.9999)

    #    for ii in range(local_lon.shape[0]):
    #        local_lon[ii,:] = np.concatenate([local_lon[ii,:][over_180] - 360.,\
    #            local_lon[ii,:][under_180]])
    #        local_lat[ii,:] = np.concatenate([local_lat[ii,:][over_180],\
    #            local_lat[ii,:][under_180]])
    #        ceres_clr_trends[ii,:] = np.concatenate([ceres_clr_trends[ii,:][over_180],\
    #            ceres_clr_trends[ii,:][under_180]])
    #        ceres_all_trends[ii,:] = np.concatenate([ceres_all_trends[ii,:][over_180],\
    #            ceres_all_trends[ii,:][under_180]])

   
    # Now average the data longitudinally
    # ----------------------------------- 
    final_naaps = np.mean(zonal_naaps, axis = 1)
    final_ceres_swclr = np.mean(zonal_ceres_swclr, axis = 1)
    final_ceres_swall = np.mean(zonal_ceres_swall, axis = 1)
    final_nsidc = np.mean(zonal_nsidc, axis = 1) 

    # Set up the figure
    plt.close('all')
    fig = plt.figure(figsize = (11,11))
    gs1 = fig.add_gridspec(nrows = 4, ncols = 3)
    ax1 = fig.add_subplot(gs1[0,0:2]) # NAAPS / Clear CERES time
    ax2 = fig.add_subplot(gs1[0,2])   # N / CC scatter
    ax3 = fig.add_subplot(gs1[1,0:2]) # NAAPS / All CERES time
    ax4 = fig.add_subplot(gs1[1,2])   # N / CA scatter
    ax5 = fig.add_subplot(gs1[2,0:2]) # NAAPS / ice time
    ax6 = fig.add_subplot(gs1[2,2])   # NAAPS / ice scatter
    ax7 = fig.add_subplot(gs1[3,0:2]) # NAAPS / ice time
    ax8 = fig.add_subplot(gs1[3,2])   # NAAPS / ice scatter

    # Plot stuff here
    # =================================
    str_dates = CERES_data['dates'][month_idx::6]
    xvals = np.array([datetime.strptime(tdate, '%Y%m') for tdate in str_dates])
    #xvals = np.arange(len(final_naaps))

   
    # Plot NAAPS versus clear 
    ax1.plot(xvals, final_naaps)
    ax1.set_ylabel('NAAPS-RA\n' + param, color = 'tab:blue')
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax12 = ax1.twinx()
    ax12.plot(xvals, final_ceres_swclr, color = 'tab:orange')
    ax12.set_ylabel('CERES Clear-sky SWF [W/m2]', color = 'tab:orange')
    ax2.scatter(final_naaps, final_ceres_swclr)
    ax2.set_xlabel('NAAPS-RA ' + param)
    ax2.set_ylabel('CERES Clear-sky SWF [W/m2]')
    ax2_result = stats.linregress(final_naaps.flatten(), final_ceres_swclr.flatten())
    corr1 = ax2_result.rvalue**2.
    ax2.set_title('r$^{2}$ = ' + str(np.round(corr1, 3)))
    #plot_figure_text(ax4, 'r$^{2}$ = ' + str(np.round(corr1, 3)), xval = None, yval = None, transform = None, \
    #    color = 'black', fontsize = 10, backgroundcolor = 'white',\
    #    halign = 'right', location = 'lower_right')
    ax1.set_title('NAAPS-RA ' + param + ' vs. CERES Clear-sky SWF', weight = 'bold')

    # Plot NAAPS versus all 
    ax3.plot(xvals, final_naaps)
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax3.set_ylabel('NAAPS-RA\n' + param, color = 'tab:blue')
    ax32 = ax3.twinx()
    ax32.plot(xvals, final_ceres_swall, color = 'tab:orange')
    ax32.set_ylabel('CERES All-sky SWF [W/m2]', color = 'tab:orange')
    ax4.scatter(final_naaps, final_ceres_swall)
    ax4.set_xlabel('NAAPS ' + param)
    ax4.set_ylabel('CERES All-sky SWF [W/m2]')
    ax4_result = stats.linregress(final_naaps.flatten(), final_ceres_swall.flatten())
    corr1 = ax4_result.rvalue**2.
    ax4.set_title('r$^{2}$ = ' + str(np.round(corr1, 3)))
    ax3.set_title('NAAPS-RA ' + param + ' vs. CERES All-sky SWF', weight = 'bold')

    # Plot NAAPS versus NSIDC
    ax5.plot(xvals, final_naaps)
    ax5.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax5.set_ylabel('NAAPS-RA\n' + param, color = 'tab:blue')
    ax52 = ax5.twinx()
    ax52.plot(xvals, final_nsidc, color = 'tab:orange')
    ax52.set_ylabel('SSMI/S Sea Ice Conc. [%]', color = 'tab:orange')
    ax6.scatter(final_naaps, final_nsidc)
    ax6.set_xlabel('NAAPS ' + param)
    ax6.set_ylabel('SSMI/S Sea Ice Conc. [%]')
    ax6_result = stats.linregress(final_naaps.flatten(), final_nsidc.flatten())
    corr1 = ax6_result.rvalue**2.
    ax6.set_title('r$^{2}$ = ' + str(np.round(corr1, 3)))
    ax5.set_title('NAAPS-RA ' + param + ' vs. SSMI/S Sea Ice Conc.', weight = 'bold')

    # Plot NAAPS versus all 
    ax7.plot(xvals, final_nsidc)
    ax7.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax7.set_ylabel('SSMI/S Sea Ice Conc. [%]', color = 'tab:blue')
    ax72 = ax7.twinx()
    ax72.plot(xvals, final_ceres_swclr, color = 'tab:orange')
    ax72.set_ylabel('CERES Clear-sky SWF [W/m2]', color = 'tab:orange')
    ax8.scatter(final_nsidc, final_ceres_swclr)
    ax8_result = stats.linregress(final_nsidc.flatten(), final_ceres_swclr.flatten())
    corr1 = ax8_result.rvalue**2.
    ax8.set_title('r$^{2}$ = ' + str(np.round(corr1, 3)))
    ax8.set_xlabel('SSMI/S Sea Ice Conc. [%]')
    ax8.set_ylabel('CERES Clear-sky SWF [W/m2]')
    ax7.set_title('SSMI/S Sea Ice Conc. vs. CERES Clear-sky SWF', weight = 'bold')
 
    plt.suptitle(str_month + ' Monthly Area Averages (2005 - 2019)\n' + \
        str(minlat) + ' $^{o}$N - ' + str(maxlat) + ' $^{o}$N\n', weight = 'bold')
 
    fig.tight_layout()

    if(save):
        outname = '_'.join(['ceres','naaps','monthly','area','avgs',param,\
            str(int(minlat)), str(int(maxlat))]) + '.png'
        fig.savefig(outname, dpi = 300)
        print('Saved image', outname)
 
    plt.show() 

    ##!#ax4_result = stats.linregress(final_naaps.flatten(), final_ceres_clr.flatten())
    ##!#ax5_result = stats.linregress(final_naaps.flatten(), final_ceres_all.flatten())
    ##!##ax8_result = stats.linregress(final_naaps.flatten(), final_ceres_clr.flatten())
    ##!#
    ##!#print(ax4_result)
    ##!#print(ax5_result)

    ##!#corr1 = ax4_result.rvalue**2.
    ##!#corr2 = ax5_result.rvalue**2.

    ##!##ax6.text(0,    80, 'r$^{2}$ = ' + str(np.round(corr1, 3)))
    ##!##ax7.text(0.05, 35, 'r$^{2}$ = ' + str(np.round(corr2, 3)), backgroundcolor = 'white')
    ##!##ax8.text(0.5,  80, 'r$^{2}$ = ' + str(np.round(corr3, 3)))

    ##!#plot_figure_text(ax4, 'r$^{2}$ = ' + str(np.round(corr1, 3)), xval = None, yval = None, transform = None, \
    ##!#    color = 'black', fontsize = 10, backgroundcolor = 'white',\
    ##!#    halign = 'right', location = 'lower_right')
    ##!#plot_figure_text(ax5, 'r$^{2}$ = ' + str(np.round(corr2, 3)), xval = None, yval = None, transform = None, \
    ##!#    color = 'black', fontsize = 10, backgroundcolor = 'white',\
    ##!#    halign = 'right', location = 'upper_right')
    ##!##plot_figure_text(ax8, 'r$^{2}$ = ' + str(np.round(corr3, 3)), xval = None, yval = None, transform = None, \
    ##!##    color = 'black', fontsize = 10, backgroundcolor = 'white',\
    ##!##    halign = 'right', location = 'upper_right')

    ##!#print("Correlation 1", ax4_result.rvalue**2.)
    ##!#print("Correlation 2", ax5_result.rvalue**2.)

    ##!#sloper = trend_type
    ##!#plot_trend_line(ax4, final_naaps.flatten(), final_ceres_clr.flatten(), \
    ##!#    color='r', linestyle = '-', slope = sloper)
    ##!#plot_trend_line(ax5, final_naaps.flatten(), final_ceres_all.flatten(), \
    ##!#    color='r', linestyle = '-', slope = sloper)
    ##!##plot_trend_line(ax8, final_naaps.flatten(), final_ceres_clr.flatten(), \
    ##!##    color='r', linestyle = '-', slope = sloper)

    ##!#ax4.set_title('NAAPS ' + param + ' Trend vs.\nCERES Clear-sky SWF Trend', fontsize = titlesize)
    ##!#ax5.set_title('NAAPS ' + param + ' Trend vs.\nCERES All-sky SWF Trend', fontsize = titlesize)
    ##!##ax8.set_title('NAAPS-RA Surface Smoke Trend vs.\nCERES Clear-sky SWF Trend', fontsize = titlesize)

    ##!#ax4.set_xlabel(param + ' Trend', fontsize = labelsize)
    ##!#ax5.set_xlabel(param + ' Trend', fontsize = labelsize)
    ##!#ax4.set_ylabel('Clear-sky SWF Trend [W m$^{-2}$]', fontsize = labelsize)
    ##!#ax5.set_ylabel('All-sky SWF Trend [W m$^{-2}$]', fontsize = labelsize)

    ##!#plot_subplot_label(ax6, '(a)', location = 'upper_left')
    ##!#plot_subplot_label(ax1, '(b)', location = 'upper_left')
    ##!#plot_subplot_label(ax2, '(c)', location = 'upper_left')
    ##!#plot_subplot_label(ax3, '(d)', location = 'upper_left')
    ##!#plot_subplot_label(ax4, '(e)', location = 'upper_left')
    ##!#plot_subplot_label(ax5, '(f)', location = 'upper_left')
    ##!##plot_subplot_label(ax6, '(f)', location = 'upper_left')
    ##!##plot_subplot_label(ax7, '(g)', location = 'upper_left')
    ##!##plot_subplot_label(ax8, '(h)', location = 'upper_left')
   
    ##!#plt.suptitle(month_obj.strftime("%B Monthly Trends (2005 - 2019)"))
 
    ##!#fig1.tight_layout()
    ##!#if(save):
    ##!#    outname = '_'.join(['ceres','naaps',param, 'monthly','trend','comp',\
    ##!#        month_obj.strftime('%b'),'combined']) + '.png'
    ##!#    fig1.savefig(outname, dpi = 300)
    ##!#    print("Saved image", outname)
    ##!#else:
    ##!#    plt.show()


def plot_compare_CERES_NSIDC_NAAPS_OMI_spatial(CERES_data, NSIDC_data, NAAPS_data, \
        OMI_data, \
        month_idx = 4, minlat = 70., trend_type = 'standard', \
        save = False):

    month_obj = datetime(year = 1,month = int(CERES_data['dates'][month_idx][4:]),day = 1)
    str_month = month_obj.strftime('%B')

    # Set up the figure
    plt.close('all')
    fig1 = plt.figure(figsize = (7,7))
    ax1  = fig1.add_subplot(2,2,1, projection = mapcrs)
    ax2  = fig1.add_subplot(2,2,2, projection = mapcrs)
    ax3  = fig1.add_subplot(2,2,3, projection = mapcrs)
    ax4  = fig1.add_subplot(2,2,4, projection = mapcrs)

    cbar_switch = True
    colorbar_label_size = 10
    ceres_trends = plotCERES_MonthTrend(CERES_data,month_idx=month_idx,save=False,\
        trend_type=trend_type,season='sunlight',minlat=minlat,return_trend=True, \
        pax = ax1, title = '')
    nsidc_trends = plotNSIDC_MonthTrend(NSIDC_data,month_idx=month_idx,save=False,\
        trend_type=trend_type,season='sunlight',minlat=minlat,\
        return_trend=True, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, \
        pax = ax2, show_pval = False, uncert_ax = None, title = '')
    naaps_trends = plotNAAPS_MonthTrend(NAAPS_data,month_idx=month_idx,\
        trend_type=trend_type,label = 'μg m$^{-3}$ per study period',\
        minlat=minlat,title = ' ', pax = ax3, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, show_pval = False, \
        uncert_ax = None, vmax = 1, return_trend = True)
    omi_trends = plotOMI_MonthTrend(OMI_data,month_idx=month_idx,\
        trend_type=trend_type,label = 'AI Trend (AI/Study Period)',\
        return_trend = True, colorbar_label_size = colorbar_label_size, \
        minlat=minlat,title = None, pax = ax4)

    ax1.set_title('Aqua CERES\nClear-Sky SWF Trends')
    ax2.set_title('SSMI/S Sea Ice\nConcentration Trends')
    ax3.set_title('NAAPS-RA Near-Surface\nSmoke Aerosol Trends')
    ax4.set_title('OMI UVAI Trends')

    plot_subplot_label(ax1, '(a)', location = 'upper_left')
    plot_subplot_label(ax2, '(b)', location = 'upper_left')
    plot_subplot_label(ax3, '(c)', location = 'upper_left')
    plot_subplot_label(ax4, '(d)', location = 'upper_left')

    plt.suptitle(month_obj.strftime('%B Trends (' + \
        CERES_data['dates'][0][:4] + ' - ' + \
        CERES_data['dates'][-1][:4] + ')'))
        
    fig1.tight_layout()
    if(save):
        outname = '_'.join(['ceres','nsidc','naaps','omi','trend','comp',\
            'spatial',month_obj.strftime('%b'), \
            CERES_data['param'].split('_')[2]]) + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()
