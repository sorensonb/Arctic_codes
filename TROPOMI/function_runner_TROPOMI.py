#!/usr/bin/env python

"""
  To access the variables in the returned dictionary:
    Latitude (2D)
      >>> trop_data['lat']

    Longitude (2D)
      >>> trop_data['lon']

    Time (1D, in datetime format)
      >>> trop_data['time']

    UVAI (2D)
      >>> trop_data['AI']

"""

from TROPOMI_Lib import *

date_str = '201908110044'
#date_str = '201807051819'
plot_TROPOMI_figure(date_str, minlat = 65., vmin = -2, vmax = 3, \
        circle_bound = True, ptitle = '', zoom = True, \
        save = False)

#filename = 'S5P_OFFL_L2__AER_AI_20190811T224359_20190812T002529_09471_01_010302_20190817T221032.nc'
#trop_data = read_TROPOMI(filename)
