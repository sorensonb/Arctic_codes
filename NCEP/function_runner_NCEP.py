#!/usr/bin/env python

"""


"""

# NOTE: data downloaded from https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html

from NCEP_Lib import *


begin_date = '20120619'
end_date = '20120623'
NCEP_data = read_NCEP(begin_date, end_date, minlat = 60.)
plot_NCEP(NCEP_data, ax = None, labelsize = 12, \
    plot_log = True, labelticksize = 10, zoom = True, vmin = None, \
    minlat = 60., circle_bound = True, vmax = None, save = False)

#date_str = '20120615'
#lat_bounds = [60, 72]
#lon_bounds = [305, 345]
#plot_NCEP_event(date_str, minlat = 60., vmin = None, vmax = None, \
#    ptitle = '', circle_bound = True, zoom = True, \
#    lat_bounds = lat_bounds, lon_bounds = lon_bounds, \
#    save = False)
