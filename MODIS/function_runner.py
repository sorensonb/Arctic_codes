#!/usr/bin/env python

"""


"""

from MODISLib import *

date_str = '202107202125'
#date_str = '202109012105'
channel1 = 1
channel2 = 5
channel3 = 31

plot_combined_figure1_v2(date_str = '202107202125', zoom = True, show_smoke = False, composite = True, \
    save=False)
#plot_true_color_satpy(date_str, ax = None, zoom = True, save = False, composite = False)
#plot_combined_figure1(date_str = date_str, zoom = True, show_smoke = True, composite = True, \
#        save=False)

#plot_figure2(save=True, composite = True, calc_radiance = True, \
#        satellite = 'modis_ch31')
#plot_scatter_OMI_CERES_figure(zoom = True, show_smoke = False, composite = True, \
#        plume_only = False, avg_pixel = True, save=True)
#plot_MODIS_detection(date_str, zoom = True, save = False)
#plot_MODIS_CERES_3panel(zoom = True, show_smoke = False, composite = True, \
#        save=False)

#compare_MODIS_3panel(date_str,channel1,channel2,channel3,zoom=True,save=True,\
#        plot_ASOS_loc = False, show_smoke = True, compare_OMI = False, \
#        compare_CERES = False, return_MODIS = False)
#plot_MODIS_temporary_4panel('202107222110', \
#    zoom = True, composite = True, show_smoke = True, save = True)
#plot_MODIS_temporary('202107222110', zoom = True, save = True)
#compare_MODIS_3panel('202107222110',31,1,5,zoom=True,save=True,\
#        plot_ASOS_loc = False, show_smoke = True)
#plot_combined_figure3(save = True)
#plot_figure2(save=True, composite = True)
#plot_figureS1(save=True, composite = True)
#plot_combined_figureS2(save=True)
