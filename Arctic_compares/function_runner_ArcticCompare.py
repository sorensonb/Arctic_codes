#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *

# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------

date_str = '201807052125'
var1 = 'NSIDC'
var2 = 'CERES_SWF'
plot_compare_scatter(date_str, var1, var2, var3 = 'MODIS_CH7', minlat = 65., \
    xmin = 1, zoom = False, save = False)
#plot_compare_colocate_spatial(date_str, minlat = 65., zoom = False, \
#    save = False)
