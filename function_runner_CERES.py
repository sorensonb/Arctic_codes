#!/usr/bin/env python

"""


"""

import gridCERESLib
from gridCERESLib import *
import sys

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

#OMI_data, CERES_data = plot_OMI_CERES_trend_compare_summer(minlat=72,\
#        ceres_type = 'sw', trend_type = 'standard', save=False)


omi_date_str = ['201708170057', \
                  '201708170236', \
                  '201708170415', \
                  '201708170554', \
                  '201708170732', \
                  '201708170911', \
                  '201708171050', \
                  '201708171229', \
                  '201708171408', \
                  '201708171547', \
                  '201708171726', \
                  '201708171905', \
                  '201708172043', \
                  '201708172222', \
                  '201708180001', \
                  '201708180140', \
                  '201708180319', \
                  '201708180458', \
                  '201708180637', \
                  '201708180816', \
                  '201708180954', \
                  '201708181133', \
                  '201708181312', \
                  '201708181451', \
                  '201708181630', \
                  '201708181809', \
                  '201708181948', \
                  '201708182127', \
                  '201708182306', \
                  '201708190044', \
                  '201708190223', \
                  '201708190402', \
                  '201708190541', \
                  '201708190720', \
                  '201708190859', \
                  '201708191038', \
                  '201708191217', \
                  '201708191355', \
                  '201708191534', \
                  '201708191713', \
                  '201708191852', \
                  '201708192031', \
                  '201708192210', \
                  '201708192349', \
                  '201908110033', \
                  '201908110212', \
                  '201908110351', \
                  '201908110530', \
                  '201908110708', \
                  '201908110847', \
                  '201908111026', \
                  '201908111205', \
                  '201908111344', \
                  '201908111523', \
                  '201908111702', \
                  '201908111841', \
                  '201908112019', \
                  '201908112158', \
                  '201908112337']

ceres_date_str = ['2017081701', \
                  '2017081702', \
                  '2017081704', \
                  '2017081705', \
                  '2017081707', \
                  '2017081709', \
                  '2017081710', \
                  '2017081713', \
                  '2017081715', \
                  '2017081716', \
                  '2017081718', \
                  '2017081720', \
                  '2017081721', \
                  '2017081723', \
                  '2017081800', \
                  '2017081802', \
                  '2017081803', \
                  '2017081804', \
                  '2017081806', \
                  '2017081808', \
                  '2017081809', \
                  '2017081811', \
                  '2017081813', \
                  '2017081814', \
                  '2017081817', \
                  '2017081819', \
                  '2017081820', \
                  '2017081822', \
                  '2017081823', \
                  '2017081900', \
                  '2017081903', \
                  '2017081904', \
                  '2017081905', \
                  '2017081907', \
                  '2017081909', \
                  '2017081910', \
                  '2017081913', \
                  '2017081914', \
                  '2017081916', \
                  '2017081918', \
                  '2017081919', \
                  '2017081920', \
                  '2017081922', \
                  '2017081923', \
                  '2019081101', \
                  '2019081103', \
                  '2019081104', \
                  '2019081105', \
                  '2019081107', \
                  '2019081108', \
                  '2019081110', \
                  '2019081112', \
                  '2019081113', \
                  '2019081115', \
                  '2019081117', \
                  '2019081118', \
                  '2019081120', \
                  '2019081121', \
                  '2019081123']

###
#OMI_date   = '201708191534'
#CERES_date = '2017081916'
#date_idx = 8
#date_idx = 9   # good
#date_idx = 10  # good
#date_idx = 11  # good
#date_idx = 12  # good
#date_idx = 13  # good
#date_idx = 24  # good
#date_idx = 25  # good
#date_idx = 26
#date_idx = 36 # kinda good
#date_idx = 38
date_idx = 39 # very very good, JZ
#date_idx = 40 # 

# 20190811
#date_idx = 44   # 2019081101
#date_idx = 45   # 2019081103
#date_idx = 46   # 2019081103

#plot_compare_OMI_CERES_hrly_grid(omi_date_str[date_idx], ceres_date_str[date_idx],  \
#date_str = '201708190000'
date_str = '201708191810'
#date_str = '201908110445'
plot_compare_OMI_CERES_hrly_grid(date_str,  \
        only_sea_ice = True, minlat = 65., skiprows = [52], \
        lat_circles = None, omi_dtype = 'JZ', zoom = True, save = False)
