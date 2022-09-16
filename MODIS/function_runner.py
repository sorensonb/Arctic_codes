#!/usr/bin/env python

"""


"""

from MODISLib import *

#date_str = '202107202125'
#date_str = '202107222110'
#date_str = '202108062025'
#date_str = '202109012105'
channel1 = 1
channel2 = 5
channel3 = 31

#date_str = '201807051950'
date_str = '201807052305'
#date_str = '201807052125'
MODIS_data = date_str
#MODIS_data = '201908110125'
write_MODIS_to_HDF5(MODIS_data, channel = 2, swath = True, \
    save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/' + date_str[:8] + '/')
write_MODIS_to_HDF5(MODIS_data, channel = 7, swath = True, \
    save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/' + date_str[:8] + '/')
sys.exit()


plot_MODIS_satpy(date_str, '7', ax = None, var = None, crs = None, \
    lons = None, lats = None, lat_lims = None, lon_lims = None, \
    vmin = None, vmax = None, ptitle = None, plabel = None, \
    labelsize = 10, colorbar = True, swath = True, zoom=False,save=False)

sys.exit()

#plot_MODIS_asos_sites(date_str, sites = None, save = False)
#plot_combined_figure1_v6(save = False)
plot_ceres_scatter(date_str, zoom=True,save=True,composite=True,\
    avg_pixel=True,plume_only=False)

sys.exit()

#plot_combined_figure1_v4(date_str = '202107222110', zoom = True, show_smoke = False, composite = True, \
#        double_fig = True, save = True)
plot_spatial_scatter(date_str, zoom = True, composite = True,\
    avg_pixel = True, plume_only = False, save = False)
sys.exit()

#plot_viewing_geometry(date_str = '202107222110', zoom = True, show_smoke = False, composite = True, \
#        save=False)
#plot_MODIS_VIIRS_SBDART(save=True, composite = True, calc_radiance = True)
#plot_spatial_scatter_wAI(date_str, zoom=True,save=True,composite=True,\
#    avg_pixel=True,plume_only=False)
#plot_combined_scatter(date_str,channel0 = 31, channel1 = 1, channel2 = 5,\
#        zoom=True,save=False,composite=True,avg_pixel=True,plume_only=False)


#plot_combined_figure1_v3(date_str = '202107202125', zoom = True, show_smoke = True, composite = True, \
#        save=False)
plot_MODIS_GOES_SBDART(save=True, composite = True, calc_radiance = True)
#plot_combined_figure1_v2(date_str = '202107202125', zoom = True, show_smoke = False, composite = True, \
#    save=True)
#plot_true_color_satpy(date_str, ax = None, zoom = True, save = False, composite = False)
#plot_combined_figure1(date_str = date_str, zoom = True, show_smoke = True, composite = True, \
#        save=True)

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
