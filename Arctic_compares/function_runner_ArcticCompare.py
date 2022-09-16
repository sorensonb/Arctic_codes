#!/usr/bin/env python

"""


"""

import Arctic_compare_lib
from Arctic_compare_lib import *

# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------

#date_str = '200804221935'
#date_str = '200804222110'
date_str = '200804222250'
#date_str = '201807051950'
#date_str = '201807052125'
#date_str = '201807052305'
#date_str = '201908110125'
#date_str = '201908110440'
var1 = 'OMI'
var2 = 'CERES_SWF'
plot_compare_OMI_CERES_MODIS_NSIDC(date_str, 7, \
    omi_dtype = 'shawn', minlat = 65., zoom = True, save = False)
sys.exit()

#plot_compare_scatter(date_str, var1, var2, var3 = 'NSIDC_LAND', minlat = 65., \
#    xmin = 1, zoom = False, save = False, trend = True)
#plot_compare_colocate_spatial(date_str, minlat = 65., zoom = False, \
#    save = False)
#cat = 'ICE_CLOUD'
#cat = 'OCEAN_CLOUD'
cat = 'LAND_CLEAR'
#plot_compare_colocate_spatial_category(date_str, cat = cat, minlat = 65., \
#    zoom = True, save = False)
trend = True

fig = plt.figure(figsize = (12,4))
ax1 = fig.add_subplot(1,3,1)
ax2 = fig.add_subplot(1,3,2)
ax3 = fig.add_subplot(1,3,3)
##ax1 = fig.add_subplot(2,3,1)
##ax2 = fig.add_subplot(2,3,4)
##ax3 = fig.add_subplot(2,3,2)
##ax4 = fig.add_subplot(2,3,5)
##ax5 = fig.add_subplot(2,3,3)
##ax6 = fig.add_subplot(2,3,6)

plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
    cat = 'ICE_CLOUD', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
    ax = ax1, colorbar = True, trend = trend, zoom = False, save = False,\
    color = 'tab:blue')
plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
    cat = 'ICE_CLEAR', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
    ax = ax1, colorbar = True, trend = trend, zoom = False, save = False,\
    color = 'tab:orange')
plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
    cat = 'OCEAN_CLOUD', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
    ax = ax2, colorbar = True, trend = trend, zoom = False, save = False, \
    color = 'tab:blue')
plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
    cat = 'OCEAN_CLEAR', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
    ax = ax2, colorbar = True, trend = trend, zoom = False, save = False, \
    color = 'tab:orange')
plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
    cat = 'LAND_CLOUD', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
    ax = ax3, colorbar = True, trend = trend, zoom = False, save = False, \
    color = 'tab:blue')
plot_compare_scatter_category(date_str, var1, var2, var3 = None, \
    cat = 'LAND_CLEAR', minlat = 65., xmin = 1, xmax = None, ymin = None, ymax = None, \
    ax = ax3, colorbar = True, trend = trend, zoom = False, save = False, \
    color = 'tab:orange')

fig.tight_layout()

outname = 'arctic_compare_scatter_6panel_' + date_str + '.png'
fig.savefig(outname, dpi = 300)
print("Saved image", outname)

plt.show()
