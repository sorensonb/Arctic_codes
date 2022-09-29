#!/usr/bin/env python

"""


"""

import NSIDCLib
from NSIDCLib import *
import sys

# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------

date_str = '20080426'
download_NSIDC_daily(date_str, save_dir = '/home/bsorenson/data/NSIDC/')
sys.exit()
##!#NSIDC_data = readNSIDC_daily(date_str, grid_data = False)
NSIDC_data = '20190811'
writeNSIDC_to_HDF5(NSIDC_data, save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/20190811/')
##!#plotNSIDC_daily_figure(date_str, minlat = 65., \
##!#    lat_circles = None, grid_data = False, zoom = False, \
##!#    vmin = None, vmax = None, circle_bound = True, \
##!#    ax = None, gridlines = False, save = False)
