"""


"""



import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime, timedelta

from siphon.simplewebservice.wyoming import WyomingUpperAir
from siphon.simplewebservice.iastate import IAStateUpperAir

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import metpy.calc as mpcalc
from metpy.interpolate import interpolate_to_grid, remove_nan_observations
from metpy.units import units
from metpy.io import add_station_lat_lon
import metpy.plots as mpplots
from metpy.plots.declarative import (BarbPlot, ContourPlot, FilledContourPlot, MapPanel,
                                     PanelContainer, PlotObs)
# Create a list of stations to pull data from
# -------------------------------------------
stations = ['UIL','SLE','MFR','REV','OAK','VBG','EDW','NKX','VEF','LKN',\
    'BOI','OTX','TFX','RIW','SLC','FGZ','TUS','EPZ','ABQ','GJT','DNR','GGW',\
    'BIS','UNR','ABR','LBF','OAX','TOP','DDC','AMA','OUN','FWD','MAF','DRT',\
    'CRP','BRO','LCH','SHV','LZK','SGF','DVN','MPX','GRB','ILX','BNA','JAN',\
    'LIX','BMX','FFC','TLH','TBW','MFL','EYW','XMR','JAX','CHS','MHX','GSO',\
    'RNK','ILN','DTX','APX','PIT','IAD','WAL','OKX','CHH','ALB','BUF','GYX',\
    'CAR','CWVK','CYZT','CYXS','CWSE','CYQD','CWPL','CYMO','CYAH','CWMW',\
    'CYZV','CYQI','CYYR']
    #'MMGM','MMCY','MMAN','MMTC','MMLP','MMMZ','MMGL',\
    #'MMZO','MMMX','MMVR','MMMD','MMUN','MYNN','MWCR','MKJP','TSDM']

contour_levels = {
    850: np.arange(1140, 1860, 30),
    700: np.arange(2640, 3360, 30),
    500: np.arange(4620, 6060, 60),
    300: np.arange(7560, 10440, 120),
    250: np.arange(8880, 11760, 120),
    200: np.arange(10560, 13440, 120)
}

def read_upper_air_iastate(date):
    try:
        data = IAStateUpperAir().request_all_data(date)
    except ValueError:
        #print("ERROR: no data available for ",date)
        #print("     Trying previous synoptic time")
        return None
        #date = date - timedelta(hours = 12)
        #data = IAStateUpperAir().request_all_data(date)
    data = add_station_lat_lon(data)
       
    return data

def read_upper_air_uwyo(date, p_lvl):

    pressures    = np.full((len(stations)),  np.nan)
    heights      = np.full((len(stations)),  np.nan)
    temperatures = np.full((len(stations)),  np.nan)
    lats         = np.full((len(stations)),  np.nan)
    lons         = np.full((len(stations)),  np.nan)
    wspd         = np.full((len(stations)),  np.nan)
    wdir         = np.full((len(stations)),  np.nan)
    
    
    for ii, stn in enumerate(stations):
    
        print(stn, p_lvl)
        ####################################################
        # Make the request (a pandas dataframe is returned).
        try:
            df = WyomingUpperAir.request_data(date, stn)
        except ValueError:
            print("No data for station",stn)
            continue
        
        ##!#####################################################
        ##!## Units are stored in a dictionary with the variable name as the key in the `units` attribute
        ##!## of the dataframe.
        ##!#print(df.units)
        ##!#
        ##!#####################################################
        ##!#print(df.units['pressure'])
        
        ##!#####################################################
        ##!## Units can then be attached to the values from the dataframe.
        ##!#pressure = df['pressure'].values * units(df.units['pressure'])
        ##!#temperature = df['temperature'].values * units(df.units['temperature'])
        ##!#dewpoint = df['dewpoint'].values * units(df.units['dewpoint'])
        ##!#u_wind = df['u_wind'].values * units(df.units['u_wind'])
        ##!#v_wind = df['v_wind'].values * units(df.units['v_wind'])
    
        if(np.min(df['pressure']) > p_lvl):
            print("ERROR: ",stn," data not high enough")
            continue
        else:
            pressures[ii]    = df['pressure'][df['pressure'] == p_lvl]
            heights[ii]      = df['height'][df['pressure'] == p_lvl]
            temperatures[ii] = df['temperature'][df['pressure'] == p_lvl]
            lats[ii]         = df['latitude'][df['pressure'] == p_lvl]
            lons[ii]         = df['longitude'][df['pressure'] == p_lvl]
            wspd[ii]         = df['speed'][df['pressure'] == p_lvl]
            wdir[ii]         = df['direction'][df['pressure'] == p_lvl]

    out_dict = {}
    out_dict['stations'] = stations
    out_dict['pressures'] = pressures
    out_dict['heights'] = heights
    out_dict['temperatures'] = temperatures
    out_dict['lats'] = lats
    out_dict['lons'] = lons
    out_dict['wspd'] = wspd
    out_dict['wdir'] = wdir

    return out_dict

def plot_upper_air_stations(data, date, p_lvl):
    obs = mpplots.PlotObs()
    obs.data = data
    obs.time = date
    obs.level = int(p_lvl) * units.hPa
    obs.fields = ['temperature','dewpoint','height']
    obs.locations = ['NW','SW','NE']
    obs.formats = [None,None,lambda v: format(v, '.0f')[:3]]
    obs.vector_field = ('u_wind','v_wind')
    obs.reduce_points = 0

    panel = mpplots.MapPanel()
    panel.layout = (1,1,1)
    panel.area = (-124,-72,20,53)
    panel.projection = 'lcc'
    panel.layers = ['coastline','borders','states','land','ocean']    
    panel.plots = [obs]

    pc = mpplots.PanelContainer()
    pc.size = (15,10)
    pc.panels = [panel]
    pc.show()

def plot_upper_air_contour(data, date, p_lvl):
    ##!#cntr2 = ContourPlot()
    ##!#cntr2.data = data.to_xarray()
    ##!#cntr2.field = 'height'
    ##!#cntr2.level = int(p_lvl) * units.hPa
    ##!#cntr2.time = date
    ##!#cntr2.contours = list(range(0, 10000, 120))
    ##!#cntr2.linecolor = 'black'
    ##!#cntr2.linestyle = 'solid'
    ##!#cntr2.clabels = True

    ##!#cfill = FilledContourPlot()
    ##!#cfill.data = data.to_xarray()
    ##!#cfill.field = 'speed'
    ##!#cfill.level = int(p_lvl) * units.hPa
    ##!#cfill.time = date
    ##!#cfill.contours = list(range(10, 201, 20))
    ##!#cfill.colormap = 'BuPu'
    ##!#cfill.colorbar = 'horizontal'
    ##!#cfill.plot_units = 'knot'

    ##!#panel = MapPanel()
    ##!#panel.area = [-125, -74, 20, 55]
    ##!#panel.projection = 'lcc'
    ##!#panel.layers = ['states', 'coastline', 'borders']
    ##!#panel.title = f'{cfill.level.m}-hPa Heights and Wind Speed at {date}'
    ##!#panel.plots = [cfill, cntr2]
    ##!#
    ##!#pc = PanelContainer()
    ##!#pc.size = (15, 15)
    ##!#pc.panels = [panel]
    ##!#
    ##!#pc.show()

    data = data.dropna()
    local_lon = data['longitude'][data['pressure'] == float(p_lvl)].to_numpy()
    local_lat = data['latitude'][data['pressure'] == float(p_lvl)].to_numpy()
    local_spd = data['speed'][data['pressure'] == float(p_lvl)].to_numpy()
    local_dir = data['direction'][data['pressure'] == float(p_lvl)].to_numpy()
    local_hgt = data['height'][data['pressure'] == float(p_lvl)].to_numpy()
    local_tmp = data['temperature'][data['pressure'] == float(p_lvl)].to_numpy()

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Objectively analyze the data
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    test_lons, test_lats, test_hgt  = remove_nan_observations(local_lon, local_lat, local_hgt)
    test_lons, test_lats, test_tmp  = remove_nan_observations(local_lon, local_lat, local_tmp)
    test_lons, test_lats, test_wspd = remove_nan_observations(local_lon, local_lat, local_spd)
    test_lons, test_lats, test_wdir = remove_nan_observations(local_lon, local_lat, local_dir)
 
    print(test_hgt)
    print(test_lons)
    print(test_lats)

    # Grid the data
    # -------------
    glon, glat, ghgt  = interpolate_to_grid(test_lons, test_lats, test_hgt, interp_type = 'linear', hres = 1)
    glon, glat, gtmp  = interpolate_to_grid(test_lons, test_lats, test_tmp, interp_type = 'linear', hres = 1)
    glon, glat, gwspd = interpolate_to_grid(test_lons, test_lats, test_wspd, interp_type = 'linear', hres = 1)
    glon, glat, gwdir = interpolate_to_grid(test_lons, test_lats, test_wdir, interp_type = 'linear', hres = 1)
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Generate the figure
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure()
    datacrs = ccrs.PlateCarree()
    mapcrs  = ccrs.AlbersEqualArea(central_longitude = -100)
    #mapcrs  = ccrs.LambertConformal()
    ax = fig.add_subplot(1,1,1, projection = mapcrs)
    
    levels = contour_levels[int(p_lvl)]
    print(levels)
    #levels = np.arange(4980, 5820, 60)
    ax.contour(glon, glat, ghgt, levels = levels, transform = datacrs)
    ax.coastlines()
    ax.add_feature(cfeature.STATES)
    ax.set_title(date.strftime('%Y-%m-%d %H:%M'))

    plt.show()
