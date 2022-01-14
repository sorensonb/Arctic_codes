#!/usr/bin/env python

"""


"""

import OMILib
from OMILib import *
import sys

sys.exit()
# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------


###minlat = 65
###plt.close('all')
###fig1 = plt.figure(figsize = (6,6))
###mapcrs = ccrs.NorthPolarStereo()
###ax0 = fig1.add_subplot(1,1,1, projection = mapcrs)
###ax0.coastlines(resolution = '50m')
###ax0.set_extent([-180,180,minlat,90],ccrs.PlateCarree())
###ax0.set_boundary(circle, transform = ax0.transAxes)
###plot_lat_circles(ax0, [80])
###plot_arctic_regions(ax0)
###ax0.gridlines()
###fig1.tight_layout()
####fig1.savefig('arctic_lat_circles_80.png',dpi=300)
###plt.show()


minlat = 65.
#OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
#OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VJZ211_2005_2020.nc', minlat = minlat)
#OMI_VSJ4   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', minlat = minlat)
#plotOMI_Compare_ClimoTrend_summer(OMI_VBS0,OMI_VJZ211,OMI_VSJ4,\
#        trend_type = 'standard', minlat=minlat,save=False)

plot_OMI_CERES_trend_compare_summer(minlat=72,\
        ceres_type = 'sw', trend_type = 'standard', save=False)

##date_str = ['201908110033',\
##            '201908110212',\
##            '201908110351',\
##            '201908110530',\
##            '201908110708',\
##            '201908110847',\
##            '201908111026',\
##            '201908111205',\
##            '201908111344',\
##            '201908111523',\
##            '201908111702',\
##            '201908111841',\
##            '201908112019',\
##            '201908112158',\
##            '201908112337',\
##            '201908120116',\
##            '201908120255',\
##            '201908120434',\
##            '201908120613',\
##            '201908120752',\
##            '201908120930',\
##            '201908121109',\
##            '201908121248',\
##            '201908121427',\
##            '201908121606',\
##            '201908121745',\
##            '201908121924',\
##            '201908122103',\
##            '201908122242']
###
#date_str = '20170818'
###    
###plotOMI_single_ground(date_str, only_sea_ice = False, minlat = 65., \
###    zoom = True, multi_panel = False, save = True)

##plotOMI_single_swath_figure(date_str[0], dtype = 'control',  \
##        only_sea_ice = False, minlat = 65., skiprows = None, lat_circles = [74], save = False)

#plot_combined_figure1(save = True)
#plot_MODIS_temporary('202107222110', zoom = True, save = True)
#plot_MODIS_temporary_4panel('202108052125', zoom = True, composite = True, show_smoke = True, save = True)
#compare_MODIS_3panel('202107222110',31,1,5,zoom=True,save=True,\
#        plot_ASOS_loc = False, show_smoke = True)
#plot_combined_figure3(save = True)
#plot_figureS1(save=True, composite = True)
#plot_combined_figureS2(save=True)
