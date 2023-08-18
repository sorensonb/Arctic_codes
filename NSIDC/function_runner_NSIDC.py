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
NSIDC_data = readNSIDC_monthly_grid_all(begin_date, end_date, \
    season, calc_month = True)

total_array    = np.full(len(NSIDC_data['dates']), np.nan)
combined_array = np.full(len(NSIDC_data['dates']), np.nan)
land_array     = np.full(len(NSIDC_data['dates']), np.nan)
ocean_array    = np.full(len(NSIDC_data['dates']), np.nan)
ice_array      = np.full(len(NSIDC_data['dates']), np.nan)
mix_array      = np.full(len(NSIDC_data['dates']), np.nan)

max_ice_for_ocean = 0
min_ice_for_ice   = 80

for tidx in range(len(NSIDC_data['dates'])):

    # Figure out the total number of grid boxes here. 
    # Land pixels also include pixels that are nonmissing in both land and ice data
    total_pixels = np.where(\
        (NSIDC_data['grid_land_conc'][tidx,:,:].mask == False) | \
        (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))[0].shape[0]
    
    combined_pixels = np.where(\
        (NSIDC_data['grid_land_conc'][tidx,:,:].mask == False) & \
        (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))
    num_combined = combined_pixels[0].shape[0]
    # The land only pixels are 
    land_only = np.where(\
        (NSIDC_data['grid_land_conc'][tidx,:,:].mask == False) & \
        (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == True))
    num_land = land_only[0].shape[0]
    oceanice_only = np.where(\
        (NSIDC_data['grid_land_conc'][tidx,:,:].mask == True) & \
        (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))
    ocean_only = np.where( \
        (((NSIDC_data['grid_land_conc'][tidx,:,:].mask == True) & \
        (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))) & 
                           (NSIDC_data['grid_ice_conc'][tidx,:,:] <= max_ice_for_ocean))
    num_ocean_only = ocean_only[0].shape[0]
    ice_only   = np.where( \
        (((NSIDC_data['grid_land_conc'][tidx,:,:].mask == True) & \
        (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))) & 
                           (NSIDC_data['grid_ice_conc'][tidx,:,:] >= min_ice_for_ice))
    num_ice_only = ice_only[0].shape[0]
    mix_only   = np.where( \
        (((NSIDC_data['grid_land_conc'][tidx,:,:].mask == True) & \
        (NSIDC_data['grid_ice_conc'][tidx,:,:].mask == False))) & 
                           ((NSIDC_data['grid_ice_conc'][tidx,:,:] > max_ice_for_ocean) &
                           (NSIDC_data['grid_ice_conc'][tidx,:,:] < min_ice_for_ice))\
        )
    num_mix_only = mix_only[0].shape[0]
  
    total_array[tidx]    = total_pixels
    combined_array[tidx] = num_combined
    land_array[tidx]     = num_land
    ocean_array[tidx]    = num_ocean_only
    ice_array[tidx]      = num_ice_only
    mix_array[tidx]      = num_mix_only

    calc_total = num_combined + num_land + num_ocean_only + num_ice_only + num_mix_only

    print(NSIDC_data['dates'][tidx], total_pixels, calc_total, num_combined, num_land, num_ocean_only,num_ice_only,num_mix_only)

pcnt_combined = ((combined_array   / total_array) * 100.)[4::6]
pcnt_land     = ((land_array       / total_array) * 100.)[4::6]
pcnt_ocean    = ((ocean_array      / total_array) * 100.)[4::6]
pcnt_ice      = ((ice_array        / total_array) * 100.)[4::6]
pcnt_mix      = ((mix_array        / total_array) * 100.)[4::6]

xvals = np.arange(NSIDC_data['data'][::6].shape[0])


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.bar(xvals, pcnt_combined, label = 'Combined')
ax.bar(xvals, pcnt_land, bottom = pcnt_combined, label = 'Land')
ax.bar(xvals, pcnt_ocean, bottom = pcnt_combined + pcnt_land, label = 'Ocean')
ax.bar(xvals, pcnt_ice, bottom = pcnt_combined + pcnt_land + pcnt_ocean, label = 'Ice')
ax.bar(xvals, pcnt_mix, bottom = pcnt_combined + pcnt_land + pcnt_ocean + pcnt_ice, \
    label = 'Mix Ice/Ocn')

ax.legend()
fig.tight_layout()

plt.show()
 
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

plotNSIDC_ClimoTrend_all(NSIDC_data,\
    trend_type = 'standard', minlat=65.,save=False)

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
