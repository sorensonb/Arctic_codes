#!/usr/bin/env python

"""


"""

from MODISLib import *

       
MODIS_data = '201908120035'
#MODIS_data = '201807052305'
swath = False
#swath = True
channel = 2 
minlat = 65.
MODIS_ch1 = read_MODIS_channel(MODIS_data, 1, swath = swath)
MODIS_ch7 = read_MODIS_channel(MODIS_data, 7, swath = swath)

cldmsk = Dataset('CLDMSK_L2_MODIS_Aqua.A2019224.0035.001.2019224173701.nc')
#cldmsk = Dataset('CLDMSK_L2_MODIS_Aqua.A2019224.0030.001.2019224173707.nc')
#cldmsk = Dataset('CLDMSK_L2_MODIS_Aqua.A2018186.2305.001.2019064032332.nc')
testmask = cldmsk['geophysical_data/Integer_Cloud_Mask'][::5,::5]
testlat  = cldmsk['geolocation_data/latitude'][::5,::5]
testlon  = cldmsk['geolocation_data/longitude'][::5,::5]
#mydmsk = Dataset('MYD35_L2.A2018186.2305.061.2018187153528.hdf')

cldmsk.close()
#cldmsk['geophysical_data/Integer_Cloud_Mask']
#

plt.close('all')
mapcrs = ccrs.NorthPolarStereo()
modis_date = MODIS_data
fig = plt.figure(figsize = (9, 9))
ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
ax4 = fig.add_subplot(2,2,4, projection = mapcrs)

zoom = False
plot_MODIS_channel(modis_date, 'true_color', swath = swath, \
    zoom = zoom, ax = ax1)
plot_MODIS_channel(modis_date, 1, swath = swath, \
    zoom = zoom, ax = ax2, vmax = 0.7)
plot_MODIS_channel(modis_date, 7, swath = swath, \
    zoom = zoom, ax = ax3, vmax = 0.4)
ax4.pcolormesh(testlon, testlat, testmask, shading = 'auto', \
    cmap = 'jet', vmin = 0, vmax = 3, transform = datacrs)
#ax1.set_extent([145, 218, 65, 80], ccrs.PlateCarree())
#ax2.set_extent([145, 218, 65, 80], ccrs.PlateCarree())
#ax3.set_extent([145, 218, 65, 80], ccrs.PlateCarree())
#ax4.set_extent([145, 218, 65, 80], ccrs.PlateCarree())
ax1.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
ax2.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
ax3.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
ax4.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
ax1.coastlines()
ax2.coastlines()
ax3.coastlines()
ax4.coastlines()
fig.tight_layout()
plt.show()
sys.exit()

write_MODIS_to_HDF5(MODIS_data, channel = 2, swath = True, \
    save_path = './', minlat = minlat, remove_empty_scans = True)

sys.exit()

date_str = '201206151425'
download_MODIS_file(date_str, dest_dir = '/home/bsorenson/data/MODIS/Aqua/')

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
