#!/usr/bin/env python

"""


"""

from MODISLib import *

fig = plt.figure(figsize = (6, 10.5))
ax1 = fig.add_subplot(4,2,1, projection = ccrs.NorthPolarStereo(central_longitude = 320))
ax2 = fig.add_subplot(4,2,2, projection = ccrs.NorthPolarStereo(central_longitude = 320))
ax3 = fig.add_subplot(4,2,3, projection = ccrs.NorthPolarStereo(central_longitude = 320))
ax4 = fig.add_subplot(4,2,4, projection = ccrs.NorthPolarStereo(central_longitude = 320))
ax5 = fig.add_subplot(4,2,5, projection = ccrs.NorthPolarStereo(central_longitude = 320))
ax6 = fig.add_subplot(4,2,6, projection = ccrs.NorthPolarStereo(central_longitude = 320))
ax7 = fig.add_subplot(4,2,7, projection = ccrs.NorthPolarStereo(central_longitude = 320))
ax8 = fig.add_subplot(4,2,8, projection = ccrs.NorthPolarStereo(central_longitude = 320))
#ax3 = fig.add_subplot(1,3,3, projection = ccrs.NorthPolarStereo(central_longitude = 320))

channel = 3 
lsize = 10
# Time 1: first column
plot_MODIS_channel('201206121525', 'true_color', swath = True, \
    zoom = True, ax = ax1, plot_borders = True, vmax = None, labelsize = lsize)
plot_MODIS_channel('201206121525', 2, swath = True, \
    zoom = True, ax = ax3, plot_borders = True, vmax = 0.8, labelsize = lsize)
plot_MODIS_channel('201206121525', 5, swath = True, \
    zoom = True, ax = ax5, plot_borders = True, vmax = None, labelsize = lsize)
plot_MODIS_channel('201206121525', 31, swath = True, \
    zoom = True, ax = ax7, plot_borders = True, vmax = None, labelsize = lsize)
# Time 2: second column
plot_MODIS_channel('201206231505', 'true_color', swath = True, \
    zoom = True, ax = ax2, plot_borders = True, vmax = None, labelsize = lsize)
plot_MODIS_channel('201206231505', 2, swath = True, \
    zoom = True, ax = ax4, plot_borders = True, vmax = 0.8, labelsize = lsize)
plot_MODIS_channel('201206231505', 5, swath = True, \
    zoom = True, ax = ax6, plot_borders = True, vmax = None, labelsize = lsize)
plot_MODIS_channel('201206231505', 31, swath = True, \
    zoom = True, ax = ax8, plot_borders = True, vmax = None, labelsize = lsize)
#plot_MODIS_channel_time_diff('201206121525','201206231505',1,zoom=True,\
#    ax = ax3, swath = True, vmin = -0.8, vmax = 0.8,\
#    circle_bound = False, plot_borders = True)
#ax1.coastlines()
#ax2.coastlines()
#ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
ax1.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
ax2.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
ax3.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
ax4.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
ax5.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
ax6.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
ax7.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
ax8.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
#ax3.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 

#plot_MODIS_satpy('201206121525', 'true_color', ax = ax1, var = None, crs = None, \
#    lons = None, lats = None, lat_lims = None, lon_lims = None, \
#    vmin = None, vmax = None, ptitle = None, plabel = None, \
#    labelsize = 10, colorbar = True, swath = True, zoom=False,save=False)
#plot_MODIS_satpy('201206231505', 'true_color', ax = ax2, var = None, crs = None, \
#    lons = None, lats = None, lat_lims = None, lon_lims = None, \
#    vmin = None, vmax = None, ptitle = None, plabel = None, \
#    labelsize = 10, colorbar = True, swath = True, zoom=False,save=False)

plt.suptitle('15:25 UTC 12-June-2012 vs 15:05 UTC 23 June 2012')

fig.tight_layout()

fig.savefig('modis_imagery_naaps_ceres_201206.png', dpi = 300)


plt.show()
#download_MODIS_file(date_str, dest_dir = '/home/bsorenson/data/MODIS/Aqua/')

sys.exit()

plot_figure2(save=False, composite = True)
#plot_combined_figure1_v6(save = True)

sys.exit()

date_str = '202107222110'
plot_ceres_scatter(date_str, zoom=True,save=True,composite=True,\
    avg_pixel=True,plume_only=False)

sys.exit()


date1 = '202107222110'
channel = 31
dt_date_str1 = datetime.strptime(date1,"%Y%m%d%H%M")
#dt_date_str2 = datetime.strptime(date2,"%Y%m%d%H%M")
filename1 = aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][dt_date_str1.strftime('%H%M')]['modis']
#filename2 = aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][dt_date_str2.strftime('%H%M')]['modis']

MODIS_data1 = read_MODIS_channel(dt_date_str1.strftime('%Y%m%d%H%M'), channel, zoom = True)

compare_data1 = nearest_grid_values(MODIS_data1)

sys.exit()

#colocate_comparison(date1, date2, channel = 31)

sys.exit()

plot_CERES_swaths(date_str = '202107222110', save = False)
sys.exit()

#date_str = '202107202125'
#date_str = '202107222110'
#date_str = '202108062025'
#date_str = '202109012105'
channel1 = 1
channel2 = 5
channel3 = 31

#date_str = '201807051950'
date_str = '201807052127'

CERES_date_str = '2008042219'
download_MODIS_swath(CERES_date_str, \
    dest_dir = '/home/bsorenson/data/MODIS/Aqua/')

#date_str = '201807052125'
MODIS_data = date_str
#MODIS_data = '201908110125'
write_MODIS_to_HDF5(MODIS_data, channel = 2, swath = True, \
    save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/' + date_str[:8] + '/')
write_MODIS_to_HDF5(MODIS_data, channel = 7, swath = True, \
    save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/' + date_str[:8] + '/')
sys.exit()


sys.exit()

#plot_MODIS_asos_sites(date_str, sites = None, save = False)

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
