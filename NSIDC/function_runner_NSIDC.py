#!/usr/bin/env python

"""


"""

import NSIDCLib
from NSIDCLib import *
import sys

# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------

##!#NSIDC_data = readNSIDC_daily('20180705')
##!#minlat = 65.
##!#writeNSIDC_to_HDF5(NSIDC_data, save_path = './', minlat = minlat, \
##!#    remove_empty_scans = True)
##!#
##!#sys.exit()

begin_date = '200504'
end_date   = '202009'
season     = 'sunlight'
minlat = 70.
NSIDC_data = readNSIDC_monthly_grid_all(begin_date, end_date, \
    season, calc_month = True, minlat = minlat, maxlat = 87.)

plotNSIDC_ClimoTrend_all(NSIDC_data,\
    trend_type = 'standard', minlat=minlat,save=False)

sys.exit()

type_dict = calc_NSIDC_sfc_types(NSIDC_data, 50, \
    use_grid_data = True)



use_area = True
def temp_plot(NSIDC_data, date_str):
    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(1,2,2, projection = mapcrs)
    plot_NSIDC_month_sfc_types(NSIDC_data, date_str, ax = ax1, save = False, \
        use_grid_data = True, max_ice_for_ocean = 0)
    plot_NSIDC_month_sfc_types(NSIDC_data, date_str, ax = ax2, save = False, \
        use_grid_data = True, max_ice_for_ocean = 20)
    
    ax1.coastlines()
    ax2.coastlines()
    ax1.set_extent([-180, 180, 70, 90], datacrs)
    ax2.set_extent([-180, 180, 70, 90], datacrs)

    ax1.set_title('Max ice for ocean = 0')
    ax2.set_title('Max ice for ocean = 20')
   
    plt.suptitle(date_str)
 
    fig.tight_layout()
    plt.show()

temp_plot(NSIDC_data, '201807')
sys.exit()


#plot_NSIDC_sfc_type_change_bar(NSIDC_data,  max_ice_for_ocean = 20, \
#    min_ice_for_ice = 80, use_area = True, use_pcnt = True, save = True)
plot_NSIDC_sfc_type_change_line(NSIDC_data, max_ice_for_ocean = 0, \
    min_ice_for_ice = 80, use_area = True, use_pcnt = True, save = False)

sys.exit() 

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
ax2 = fig.add_subplot(1,2,2, projection = mapcrs)

ax1.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], NSIDC_data['grid_land_conc'][tidx,:,:], \
    transform = datacrs, shading = 'auto')
ax2.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], NSIDC_data['grid_ice_conc'][tidx,:,:], \
    transform = datacrs, shading = 'auto')

ax1.coastlines()
ax2.coastlines()

ax1.set_extent([-180, 180, 65, 90], datacrs)
ax2.set_extent([-180, 180, 65, 90], datacrs)

plt.suptitle(NSIDC_data['dates'][tidx])

fig.tight_layout()

plt.show()

sys.exit()
plotNSIDC_MonthTrend(NSIDC_data,month_idx=4,save=False,\
    trend_type='theil-sen',season='',minlat=70.,return_trend=False, \
    colorbar = True, colorbar_label_size = None,title = None, \
    pax = None, show_pval = False, uncert_ax = None)

sys.exit()


fig = plt.figure(figsize = (10, 11))
ax1 = fig.add_subplot(3,3,1, projection = mapcrs)
ax2 = fig.add_subplot(3,3,2, projection = mapcrs)
ax3 = fig.add_subplot(3,3,3, projection = mapcrs)
ax4 = fig.add_subplot(3,3,4, projection = mapcrs)
ax5 = fig.add_subplot(3,3,5, projection = mapcrs)
ax6 = fig.add_subplot(3,3,6, projection = mapcrs)
ax7 = fig.add_subplot(3,3,7, projection = mapcrs)
ax8 = fig.add_subplot(3,3,8, projection = mapcrs)
ax9 = fig.add_subplot(3,3,9, projection = mapcrs)

ax1.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    NSIDC_data['MONTH_CLIMO'][0,:,:], shading = 'auto', \
    transform = datacrs)
ax1.coastlines()
ax1.set_extent([-180,180,65,90], datacrs)

ax2.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    NSIDC_data['MONTH_CLIMO'][1,:,:], shading = 'auto', \
    transform = datacrs)
ax2.coastlines()
ax2.set_extent([-180,180,65,90], datacrs)

ax3.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    NSIDC_data['MONTH_CLIMO'][2,:,:], shading = 'auto', \
    transform = datacrs)
ax3.coastlines()
ax3.set_extent([-180,180,65,90], datacrs)

ax4.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    NSIDC_data['MONTH_CLIMO'][3,:,:], shading = 'auto', \
    transform = datacrs)
ax4.coastlines()
ax4.set_extent([-180,180,65,90], datacrs)

ax5.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    NSIDC_data['MONTH_CLIMO'][4,:,:], shading = 'auto', \
    transform = datacrs)
ax5.coastlines()
ax5.set_extent([-180,180,65,90], datacrs)

ax6.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    NSIDC_data['MONTH_CLIMO'][5,:,:], shading = 'auto', \
    transform = datacrs)
ax6.coastlines()
ax6.set_extent([-180,180,65,90], datacrs)

ax7.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    NSIDC_data['grid_ice_conc'][34,:,:], shading = 'auto', \
    transform = datacrs)
ax7.coastlines()
ax7.set_extent([-180,180,65,90], datacrs)

ax8.pcolormesh(NSIDC_data['grid_lon'], NSIDC_data['grid_lat'], \
    NSIDC_data['grid_ice_conc'][12,:,:], shading = 'auto', \
    transform = datacrs)
ax8.coastlines()
ax8.set_extent([-180,180,65,90], datacrs)

ax9.pcolormesh(NSIDC_data['lon'], NSIDC_data['lat'], \
    NSIDC_data['data'][12,:,:], shading = 'auto', \
    transform = datacrs)
ax9.coastlines()
ax9.set_extent([-180,180,65,90], datacrs)

fig.tight_layout()
plt.show()

sys.exit()

local_date = datetime(2005,4,1)
download_NSIDC_monthly(local_date.strftime('%Y%m'))

sys.exit()

begin_date = datetime(2005,4,1)
end_date   = datetime(2020,9,1)

local_date = begin_date
while(local_date <= end_date):
    print(local_date)
    local_date = local_date + relativedelta(months = 1)

    download_NSIDC_monthly(local_date.strftime('%Y%m'))

sys.exit()

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
