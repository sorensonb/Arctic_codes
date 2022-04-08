#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>  - 2020/02/19
      Written (based on Wyoming_request.py and Skew-T_Layout.py)

"""

# Copyright (c) 2017 Siphon Contributors.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
"""
Wyoming Upper Air Data Request
==============================

This example shows how to use siphon's `simplewebswervice` support to create a query to
the Wyoming upper air archive.
"""

import sys
import matplotlib.pyplot as plt
from datetime import datetime
from siphon.simplewebservice.wyoming import WyomingUpperAir

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import metpy.calc as mpcalc
#from metpy.plots import add_metpy_logo, interpolate_to_grid, remove_nan_observations
from metpy.interpolate import interpolate_to_grid, remove_nan_observations
from metpy.units import units

from met_data_lib import *

if(len(sys.argv)!=2):
    print("SYNTAX: python upper_air_plotter.py pressure")
    print("        pressure: formatted as 925, 850, 700, 500, etc.")
    sys.exit()

p_lvl = float(sys.argv[1])

####################################################
# Create a datetime object for the sounding and string of the station identifier.
#date = datetime(2017, 9, 10, 12)
#station = 'BIS'
date = datetime.utcnow().replace(microsecond=0,second=0,minute=0)
if(date.hour>=12):
    print("retrieving for ",date)
    date = datetime.utcnow().replace(microsecond=0,second=0,minute=0,hour=12)
else:
    print("retrieving for ",date)
    date = datetime.utcnow().replace(microsecond=0,second=0,minute=0,hour=0)


data = read_upper_air_iastate(date)
#if(data == None):
#    print("ERROR: no data available for ",date)
#    print("     Trying previous synoptic time")
#    date = date - timedelta(hours = 12)
#    data = read_upper_air_iastate(date)
#plot_upper_air_stations(data, date, p_lvl)
plot_upper_air_contour(data, date, p_lvl)

##!#data_dict = read_upper_air_uwyo(date, p_lvl)
##!#
##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
##!##
##!## Objectively analyze the data
##!##
##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
##!#test_lons, test_lats, test_hgt  = remove_nan_observations(data_dict['lons'], data_dict['lats'], data_dict['heights'])
##!#test_lons, test_lats, test_tmp  = remove_nan_observations(data_dict['lons'], data_dict['lats'], data_dict['temperatures'])
##!#test_lons, test_lats, test_wspd = remove_nan_observations(data_dict['lons'], data_dict['lats'], data_dict['wspd'])
##!#test_lons, test_lats, test_wdir = remove_nan_observations(data_dict['lons'], data_dict['lats'], data_dict['wdir'])
##!#
##!## Grid the data
##!## -------------
##!#glon, glat, ghgt  = interpolate_to_grid(test_lons, test_lats, test_hgt, interp_type = 'linear', hres = 1)
##!#glon, glat, gtmp  = interpolate_to_grid(test_lons, test_lats, test_tmp, interp_type = 'linear', hres = 1)
##!#glon, glat, gwspd = interpolate_to_grid(test_lons, test_lats, test_wspd, interp_type = 'linear', hres = 1)
##!#glon, glat, gwdir = interpolate_to_grid(test_lons, test_lats, test_wdir, interp_type = 'linear', hres = 1)
##!#
##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
##!##
##!## Generate the figure
##!##
##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
##!#
##!## Create a new figure. The dimensions here give a good aspect ratio
##!#fig = plt.figure()
##!#datacrs = ccrs.PlateCarree()
##!#mapcrs  = ccrs.AlbersEqualArea(central_longitude = -100)
##!##mapcrs  = ccrs.LambertConformal()
##!#ax = fig.add_subplot(1,1,1, projection = mapcrs)
##!#
##!#levels = np.arange(4980, 5820, 60)
##!#ax.contour(glon, glat, ghgt, levels = levels, transform = datacrs)
##!#ax.coastlines()
##!#ax.add_feature(cfeature.STATES)
##!#ax.set_title(date.strftime('%Y-%m-%d %H:%M'))
##!#
##!## Show the plot
##!#plt.show()
