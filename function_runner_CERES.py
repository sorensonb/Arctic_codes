#!/usr/bin/env python

"""


"""

import gridCERESLib
from gridCERESLib import *
import sys

data_dt = '2018070409'
fig = plt.figure()
minlat = 65.
plotCERES_hrly_figure(data_dt, 'SWF',  \
    minlat = minlat, lat_circles = None, ax = None, title = 'SWF',\
    grid_data = True, zoom = False, vmax = 450, vmin = None, save = False)
sys.exit()


def local_test_func(infile, data_dt, minlat, maxlat, minlon, maxlon, sizer = 120, \
        vmin1 = 150, vmax1 = 250, vmin2 = 300, vmax2 = 370):

    str_fmt = "%Y%m%d%H%M"
    minute_adder = 30
    minute_subtracter = 30

    dt_data_begin = datetime.strptime(data_dt,str_fmt) \
        - relativedelta(minutes = minute_subtracter)
    dt_data_end   = datetime.strptime(data_dt,str_fmt) \
        + relativedelta(minutes = minute_adder)

    base_date = datetime(year=1970,month=1,day=1)

    data = Dataset(infile,'r')
    time  = data.variables['time'][:]
    if('NOAA20' in infile):
        lat = data['Time_and_Position/instrument_fov_latitude'][:]
        lon = data['Time_and_Position/instrument_fov_longitude'][:]
        sza = data['Viewing_Angles/solar_zenith_angle'][:]
        vza = data['Viewing_Angles/view_zenith_angle'][:]
        azm = data['Viewing_Angles/relative_azimuth_angle'][:]
        lon[lon>179.99] = -360.+lon[lon>179.99]
        lwf = data['TOA_and_Surface_Fluxes/toa_longwave_flux'][:] 
        swf = data['TOA_and_Surface_Fluxes/toa_shortwave_flux'][:] 
    else:
        lat   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        lon   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        sza   = data.variables['CERES_solar_zenith_at_surface'][:]
        vza   = data.variables['CERES_viewing_zenith_at_surface'][:]
        azm   = data.variables['CERES_relative_azimuth_at_surface'][:]
        lon[lon>179.99] = -360.+lon[lon>179.99]
        lwf  = data.variables['CERES_LW_TOA_flux___upwards'][:]
        swf  = data.variables['CERES_SW_TOA_flux___upwards'][:]

    data.close()
    total_times = np.array([base_date + relativedelta(days = ttime) \
        for ttime in time])

    test_time = total_times[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_ftim = time[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_swf  = swf[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lwf  = lwf[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lat  = lat[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_lon  = lon[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_sza  = sza[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_vza  = vza[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]
    test_azm  = azm[np.where((total_times >= dt_data_begin) & (total_times <= dt_data_end))]

    keep_idxs = np.where( (test_lat >= minlat) & (test_lat <= maxlat) & \
            (test_lon >= minlon) & (test_lon <= maxlon))

    final_lats = test_lat[keep_idxs]
    final_lons = test_lon[keep_idxs]
    final_swfs = test_swf[keep_idxs]
    final_lwfs = test_lwf[keep_idxs]

    plt.close('all')
    fig = plt.figure(figsize = (9, 5))
    ax1 = fig.add_subplot(1,2,1, projection = ccrs.LambertConformal(\
        central_longitude = (minlon + maxlon) / 2))
    ax2 = fig.add_subplot(1,2,2, projection = ccrs.LambertConformal(\
        central_longitude = (minlon + maxlon) / 2))

    scat = ax1.scatter(final_lons, final_lats, c = final_swfs, transform = datacrs, \
        cmap = 'plasma', s = sizer, vmin = vmin1, vmax = vmax1)
    cbar = fig.colorbar(scat, ax = ax1)
    ax1.add_feature(cfeature.STATES)
    ax1.coastlines()
    ax1.set_extent([minlon, maxlon, minlat, maxlat], datacrs)

    scat = ax2.scatter(final_lons, final_lats, c = final_lwfs, transform = datacrs, \
        cmap = 'plasma', s = sizer, vmin = vmin2, vmax = vmax2)
    cbar = fig.colorbar(scat, ax = ax2)
    ax2.add_feature(cfeature.STATES)
    ax2.coastlines()
    ax2.set_extent([minlon, maxlon, minlat, maxlat], datacrs)

    plt.show()
    
#date_str = '2017081421'
#date_str = '2019081023'
####date_str = '2021072309'
####minlat = 30.0
###param = 'SWF'
###date_str = '2021072321'
#CERES_data_hrly = readgridCERES_hrly_grid(date_str[:10], param, \
#    satellite = 'NOAA20', minlat = 30.)
###sys.exit()
#infile = '/home/bsorenson/data/CERES/SSF_Level2/Aqua/modis_comp/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072208-2021072222.nc'
#infile = '/home/bsorenson/data/CERES/SSF_Level2/Terra/CERES_SSF_Terra-XTRK_Edition4A_Subset_2021072303-2021072321.nc'
infile = '/home/bsorenson/data/CERES/SSF_Level2/NOAA20/CERES_SSF_NOAA20-XTRK_Edition1B_Subset_2021072307-2021072323.nc'
#infile = '/home/bsorenson/data/CERES/SSF_Level2/Aqua/modis_comp/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072309-2021072321.nc'

#date_str = '2021072321'
#vmin1 = 170
#vmax1 = 240
#vmin2 = 300
#vmax2 = 370
#sizer = 900

#date_str = '2021072318'
#date_str = '2021072309'   # Aqua
#date_str = '2021072306'    # Terra
date_str = '202107231040'   # NOAA20
vmin1 = 150
vmax1 = 250
vmin2 = 270
vmax2 = 310
sizer = 50
dlat = 5.

#local_test_func(infile, date_str, 39.5 - dlat, 42. + dlat, -122. - dlat, -119.5 + dlat, sizer = sizer, \
#    vmin1 = vmin1, vmax1 = vmax1, vmin2 = vmin2, vmax2 = vmax2)
"""
tester = readgridCERES_hrly_grid(date_str[:10], 'SWF', satellite = 'NOAA20', minlat = 30)
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection = ccrs.LambertConformal(central_longitude = -121))
mesh = ax.pcolormesh(tester['lon'], tester['lat'], tester['lwf'], shading = 'auto', transform = datacrs, cmap = 'plasma', vmin = 270)
cbar = fig.colorbar(mesh, ax = ax, label = 'TOA LWF [W/m2]')
ax.set_title('NOAA-20 CERES TOA Flux\n23 July 2021 10:00 UTC')
ax.coastlines()
ax.add_feature(cfeature.STATES)
ax.set_extent([-122,-119.5,39.5,42], datacrs)
fig.tight_layout()
fig.savefig('ceres_noaa20_lwf_2021072310.png', dpi = 200)
plt.show()
"""
#plotCERES_hrly_figure(date_str, 'SWF',  \
#    satellite = 'Terra',
#    only_sea_ice = False, minlat = 30., \
#    lat_circles = None, grid_data = True, zoom = False, \
#    vmax = None, vmin = None, save = False)
#    #vmax = 450, vmin = None, save = False)
sys.exit()


start_date = '200504'
end_date   = '201909'
season     = 'sunlight'

# NOTE: gridded CERES SSF_1Deg data are the same along
#       latitude lines, and along wider ranges farther
#       north. The result is "pixels" in the data that 
#       appear roughly the same size. 
minlat = 60.
CERES_swclr = readgridCERES(start_date,end_date,'toa_sw_clr_mon',\
    satellite = 'Aqua',minlat=minlat,calc_month = True,season = 'sunlight')
CERES_swall = readgridCERES(start_date,end_date,'toa_sw_all_mon',\
    satellite = 'Aqua',minlat=minlat,calc_month = True,season = 'sunlight')

NSIDC_data = readNSIDC_monthly_grid_all(start_date, end_date, \
    season, minlat = minlat, calc_month = True)

NSIDC_data['grid_lat'] += 0.5

param = 'smoke'
NAAPS_data = readgridNAAPS_monthly(start_date,end_date,param,minlat=minlat,\
             calc_month = True,season = 'sunlight', calc_yearly_accumulation = True)


param = 'smoke_conc_sfc_yearly'
#plot_compare_CERES_NAAPS_monthly(CERES_swclr, CERES_swall, \
#    NSIDC_data, NAAPS_data, param = param, \
#    month_idx = 4, minlat = minlat, trend_type = 'standard', \
#    save = False)


minlat = 60.5
maxlat = 70.5

plot_compare_CERES_NAAPS_monthly_timeseries(CERES_swclr, CERES_swall, \
    NSIDC_data, NAAPS_data, param = param, \
    month_idx = 4, minlat = minlat, maxlat = maxlat, trend_type = 'standard', \
    save = False)

##!#local_CERES_swclr = np.ma.masked_where(CERES_swclr['data'] == -999., CERES_swclr['data'])
##!#local_CERES_swall = np.ma.masked_where(CERES_swall['data'] == -999., CERES_swall['data'])
##!#
##!#naaps_idx = np.where((NAAPS_data['lats'][:,0] >= minlat) & (NAAPS_data['lats'][:,0] <= maxlat))
##!#ceres_idx = np.where((CERES_swclr['lat'][:,0] >= minlat) & (CERES_swclr['lat'][:,0] <= maxlat))
##!#nsidc_idx = np.where((NSIDC_data['grid_lat'][:,0] >= minlat) & (NSIDC_data['grid_lat'][:,0] <= maxlat))
##!#
##!#zonal_naaps = np.mean(NAAPS_data['smoke_totsink_yearly'][4::6,naaps_idx[0],:], axis = 2)
##!#zonal_ceres_swclr = np.mean(local_CERES_swclr[4::6,naaps_idx[0],:], axis = 2)
##!#zonal_ceres_swall = np.mean(local_CERES_swall[4::6,naaps_idx[0],:], axis = 2)
##!#zonal_nsidc = np.mean(NSIDC_data['grid_ice_conc'][4::6,naaps_idx[0],:], axis = 2)
##!#
##!#final_naaps = np.mean(zonal_naaps, axis = 1)
##!#final_ceres_swclr = np.mean(zonal_ceres_swclr, axis = 1)
##!#final_ceres_swall = np.mean(zonal_ceres_swall, axis = 1)
##!#final_nsidc = np.mean(zonal_nsidc, axis = 1)

sys.exit()


# ----------------------------------------------------
# Set up the overall figure
# ----------------------------------------------------

date_str = '2018070523'
minlat = 65.
CERES_grid_hrly =  readgridCERES_hrly_grid(date_str, 'swf', \
    minlat=minlat)

write_CERES_hrly_grid_to_HDF5(CERES_grid_hrly, save_path = './', \
    minlat = 65., remove_empty_scans = True)

sys.exit()

#plotCERES_MonthTrend_AllClr('200504', '202009', 4, \
#    minlat = 70, satellite = 'Aqua', save = True)
start_date = '200504'
end_date   = '202009'
season     = 'sunlight'

# NOTE: gridded CERES SSF_1Deg data are the same along
#       latitude lines, and along wider ranges farther
#       north. The result is "pixels" in the data that 
#       appear roughly the same size. 
minlat = 70.
CERES_swclr = readgridCERES(start_date,end_date,'toa_sw_clr_mon',\
    satellite = 'Aqua',minlat=minlat,calc_month = True,season = 'sunlight')
CERES_swall = readgridCERES(start_date,end_date,'toa_sw_all_mon',\
    satellite = 'Aqua',minlat=minlat,calc_month = True,season = 'sunlight')

NSIDC_data = readNSIDC_monthly_grid_all(start_date, end_date, \
    season, minlat = minlat, calc_month = True)

NAAPS_data = readgridNAAPS_NCDF(infile=home_dir + \
    '/Research/NAAPS/naaps_grid_smoke_conc_sfc_2005_2020_noAI.nc',\
    start_date = 200504, end_date = int(end_date), calc_month = True, \
    minlat = minlat)

#OMI_VJZ211 = readOMI_NCDF(infile = '/home/bsorenson/Research/OMI/omi_ai_VJZ211_2005_2020.nc', minlat = minlat)
OMI_data   = readOMI_NCDF(infile = \
    '/home/bsorenson/Research/OMI/omi_ai_VSJ4_2005_2020.nc', \
    minlat = minlat, end_date = int(end_date))

#plot_compare_CERES_NSIDC_NAAPS_OMI(CERES_swall, NSIDC_data, NAAPS_data, \
#    OMI_data, month_idx = 3, minlat = minlat, trend_type = 'theil-sen', \
#    save = True)
#
#sys.exit()

for ii in range(5):
    #plot_compare_CERES_NSIDC_NAAPS_OMI_spatial(CERES_swclr, NSIDC_data, NAAPS_data, \
    #    OMI_data, month_idx = ii, minlat = minlat, trend_type = 'theil-sen', \
    #    save = True)
    plot_compare_CERES_NSIDC_NAAPS_OMI(CERES_swall, NSIDC_data, NAAPS_data, \
        OMI_data, month_idx = ii, minlat = minlat, trend_type = 'theil-sen', \
        save = True)

sys.exit()

plotCERES_ClimoTrend_all(CERES_swclr,\
    trend_type = 'standard', minlat=70.,save=True)

sys.exit()

start_date = '200503'
end_date   = '202009'
CERES_swall = readgridCERES(start_date,end_date,'toa_sw_all_mon',\
    satellite = 'Aqua',minlat=60.5,calc_month = True,season = 'sunlight')
CERES_lwall = readgridCERES(start_date,end_date,'toa_lw_all_mon',\
    satellite = 'Aqua',minlat=60.5,calc_month = True,season = 'sunlight')
CERES_swclr = readgridCERES(start_date,end_date,'toa_sw_clr_mon',\
    satellite = 'Aqua',minlat=60.5,calc_month = True,season = 'sunlight')
CERES_lwclr = readgridCERES(start_date,end_date,'toa_lw_clr_mon',\
    satellite = 'Aqua',minlat=60.5,calc_month = True,season = 'sunlight')
#CERES_all  = readgridCERES(start_date,end_date,'toa_sw_all_mon',\
#    satellite = 'Aqua',minlat=60.5,calc_month = True,season = 'sunlight')

fig1 = plt.figure(figsize = (10,10))
ax1 = fig1.add_subplot(2,2,1, projection = mapcrs)
ax2 = fig1.add_subplot(2,2,2, projection = mapcrs)
ax3 = fig1.add_subplot(2,2,3, projection = mapcrs)
ax4 = fig1.add_subplot(2,2,4, projection = mapcrs)

plotCERES_MonthTrend(CERES_swclr,month_idx=4,save=False,\
    trend_type='standard',season='sunlight',minlat=65.,return_trend=False, \
    pax = ax1)
plotCERES_MonthTrend(CERES_lwclr,month_idx=4,save=False,\
    trend_type='standard',season='sunlight',minlat=65.,return_trend=False, \
    pax = ax2)
plotCERES_MonthTrend(CERES_swall,month_idx=4,save=False,\
    trend_type='standard',season='sunlight',minlat=65.,return_trend=False, \
    pax = ax3)
plotCERES_MonthTrend(CERES_lwall,month_idx=4,save=False,\
    trend_type='standard',season='sunlight',minlat=65.,return_trend=False, \
    pax = ax4)

fig1.tight_layout()
plt.show()

sys.exit()

write_CERES_L2_to_HDF5('20170816', 'Aqua', save_path = './')

sys.exit()

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


#date_str = '2008042219' # GOOD
#date_str = '2008042221' # GOOD
#date_str = '2008042222' # GOOD
#date_str = '2016051520' # MEDIOCRE
#date_str = '2016051522' # MEDIOCRE
#date_str = '2016051523' # MEDIOCRE
#date_str = '2016051622' # MEDIOCRE
#date_str = '2018070519' # GOOD
#date_str = '2018070521' # GOOD
#date_str = '2018070523' # GOOD
#date_str = '2019081022' # GOOD
#date_str = '2019081023' # GOOD
#date_str = '2019081101' # GOOD
#date_str = '2019081104' # GOOD
#date_str = '2006072601' # GOOD
#date_str = '2006072523' # GOOD
#date_str = '2006072602' # GOOD
#date_str = '2006072604' # GOOD
#date_str = '2006072606' # GOOD
#date_str = '2017081616' # GOOD
#date_str = '2017081617' # GOOD
#date_str = '2017081619' # GOOD
#date_str = '2017081715' # GOOD
#date_str = '2017081716' # GOOD
#date_str = '2017081718' # GOOD
#date_str = '2017081720' # GOOD
#date_str = '2017081721' # GOOD
#date_str = '2017081814' # GOOD
#date_str = '2017081815' # GOOD
#date_str = '2017081817' # GOOD
#date_str = '2017081819' # GOOD
#date_str = '2017081820' # GOOD
#date_str = '2017081914' # GOOD
#date_str = '2017081916' # GOOD
#date_str = '2017081918' # GOOD

#date_str = '2006072701'
#plotCERES_hrly_figure(date_str, 'SWF',  \
#    only_sea_ice = False, minlat = 65., \
#    lat_circles = None, grid_data = True, zoom = False, \
#    vmax = 450, vmin = None, save = False)
#sys.exit()

CERES_data = '20170820'
end_str    = '20170831'
pvar = 'alb_clr'
#plotCERES_daily(CERES_data, pvar, end_str = end_str, satellite = 'Aqua',  \
#    only_sea_ice = False, minlat = 65., avg_data = True, \
#    lat_circles = None, ax = None, save = False, \
#    circle_bound = True, colorbar = True)

plotCERES_daily_allsat('20170801', pvar, end_str = '20170815', \
    only_sea_ice = False, minlat = 65., avg_data = True, \
    lat_circles = None, save = True, \
    circle_bound = True, colorbar = True)
#CERES_data1 = readgridCERES_daily('20170801',end_str = '20170815', satellite = 'Aqua',minlat=70.5)
#CERES_data2 = readgridCERES_daily('20170801',end_str = '20170815', satellite = 'Terra',minlat=70.5)
#CERES_data3 = readgridCERES_daily('20170801',end_str = '20170815', satellite = 'SuomiNPP',minlat=70.5)
#CERES_data2 = readgridCERES_daily('20170820',end_str = '20170831', satellite = 'Aqua',minlat=70.5)
sys.exit()

param = 'SWF'
CERES_data_hrly = readgridCERES_hrly_grid(date_str, param, \
    satellite = 'Aqua', minlat = minlat)


CERES_grid_hrly = date_str
#CERES_grid_hrly = '2019081101'
write_CERES_hrly_grid_to_HDF5(CERES_grid_hrly, \
    save_path = '/home/bsorenson/Research/Arctic_compares/comp_data/20180705/')
sys.exit()



# 20190811
#date_idx = 44   # 2019081101
#date_idx = 45   # 2019081103
#date_idx = 46   # 2019081103

#plot_compare_OMI_CERES_hrly_grid(omi_date_str[date_idx], ceres_date_str[date_idx],  \
OMI_date   = '201708191713'
CERES_date = '2017081918'
#date_str = '201708190000'
date_str = '201708191810'
#date_str = '201908110445'
#plot_compare_OMI_CERES_hrly_grid(date_str,  \
#        only_sea_ice = False, minlat = 65., skiprows = [52], \
#        no_ice = False, lat_circles = None, omi_dtype = 'control', zoom = True, save = False)
plot_compare_OMI_CERES_hrly_grid_2case('201908110445','201708191810', minlat=65,max_AI = -200.,\
        omi_dtype = 'control', ceres_dtype = 'swf', only_sea_ice = False, \
        only_ice = False, no_ice = False,  skiprows = [52], save=True, \
        composite = False, zoom = True, lat_circles = None, \
        show_scatter = False)

