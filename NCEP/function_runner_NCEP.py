#!/usr/bin/env python

"""


"""

# NOTE: data downloaded from https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html

from NCEP_Lib import *


#begin_date = '20170801'
#end_date = '20170815'
#NCEP_data = read_NCEP(begin_date, end_date, minlat = 65.)

date_str = '20170816'
lat_bounds = [78, 85]
lon_bounds = [180, 340]
plot_NCEP_event(date_str, minlat = 65., vmin = None, vmax = None, \
    ptitle = '', circle_bound = True, zoom = True, \
    lat_bounds = lat_bounds, lon_bounds = lon_bounds, \
    save = False)
