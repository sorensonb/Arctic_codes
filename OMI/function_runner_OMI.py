#!/usr/bin/env python

"""


"""

import OMILib
from OMILib import *

minlat = 65.
OMI_VBS0   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VBS0_2005_2020.nc', minlat = minlat)
OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VJZ211_2005_2020.nc', minlat = minlat)
OMI_VSJ4   = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', minlat = minlat)

plotOMI_Compare_ClimoTrend_summer(OMI_VBS0,OMI_VJZ211,OMI_VSJ4,\
        trend_type = 'standard', minlat=minlat,save=True)

#plot_OMI_CERES_trend_compare_summer(minlat=72,\
#        ceres_type = 'lw', trend_type = 'theil-sen', save=False)

#plot_combined_figure1(save = True)
#plot_MODIS_temporary('202107222110', zoom = True, save = True)
#plot_MODIS_temporary_4panel('202108052125', zoom = True, composite = True, show_smoke = True, save = True)
#compare_MODIS_3panel('202107222110',31,1,5,zoom=True,save=True,\
#        plot_ASOS_loc = False, show_smoke = True)
#plot_combined_figure3(save = True)
#plot_figureS1(save=True, composite = True)
#plot_combined_figureS2(save=True)
