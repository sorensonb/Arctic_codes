#!/usr/bin/env python

"""


"""

from NAAPSLib import *

#x = 1010
#
#for ii in range(5):
#    str_x = sorted(str(x))[::-1]
#    str_y = str_x[::-1]
#    if(len(str_y) < 4):
#        filler = ['0'] * (4 - len(str_y))
#        str_y = str_y + filler
#    x = int(''.join(str_x))
#    y = int(''.join(str_y))
#    z = x - y
#    print(x,'   ',y,'   ',z)
#
#    x = z
#
#sys.exit()

save_flag = True
minlat = 65.
date_str = '20080422'
shawn_path = home_dir + '/data/OMI/shawn_files/ltc3_old/'
#var = 'smoke_conc_sfc'
#var = 'smoke_drysink'
#var = 'smoke_wetsink'

shawn_files = glob(shawn_path + '200804*')
dt_dates = np.array([datetime.strptime(sfile.strip().split('/')[-1], \
    '%Y%m%d%H%M') for sfile in shawn_files])
date_strs = np.array([dtd.strftime('%Y%m%d%H%M') for dtd in dt_dates])

#sys.exit()

base_date = datetime(2008,4,20,18)
for ii in range(16):
    base_date = base_date + timedelta(hours = 6)
    date_str = base_date.strftime('%Y%m%d%H')    

    # Load and plot Shawn averaged data
    # Figure out which swaths are within the range
    # --------------------------------------------
    min_window = base_date - timedelta(hours = 3)
    max_window = base_date + timedelta(hours = 3)

    in_times = date_strs[(dt_dates >= min_window) & \
        (dt_dates <= max_window)]


    OMI_data = readOMI_swath_shawn_old(in_times[0], latmin = minlat, \
        resolution = 1.00, shawn_path = shawn_path)
    local_AI = np.full((len(in_times), OMI_data['AI'].shape[0], \
        OMI_data['AI'].shape[1]), np.nan)
    local_AI[0,:,:] = OMI_data['AI'][:,:]

    for ii in range(1, len(in_times)):
        OMI_local = readOMI_swath_shawn_old(in_times[ii], latmin = minlat,
            resolution = 1.00, shawn_path = shawn_path)
        local_AI[ii,:,:] = OMI_local['AI'][:,:]

    mask_AI = np.ma.masked_where(local_AI == 0., local_AI)
    mask_AI = np.ma.masked_invalid(mask_AI)
    mask_AI = np.nanmean(mask_AI, axis = 0)

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = mapcrs)
    mesh = ax.pcolormesh(OMI_data['LON'], OMI_data['LAT'], mask_AI.T, \
        transform = datacrs, shading = 'auto', cmap = 'jet', \
        vmin = -2, vmax = 3)
    cbar = plt.colorbar(mesh, ax = ax, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label('UVAI', size = None, weight = None)
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_extent([-180,180,minlat, 90], datacrs)
    ax.coastlines()
    ax.set_title('OMI Perturbed UVAI\n' + \
        min_window.strftime('%Y-%m-%d %H UTC') + ' - '  + \
        max_window.strftime('%Y-%m-%d %H UTC'))
    fig.tight_layout()
    outname = 'omi_assim_window_' + min_window.strftime('%Y%m%d%H') + \
        '_' + max_window.strftime('%Y%m%d%H') + '.png'
    fig.savefig(outname, dpi = 300)
    print("Saved image", outname)

    plot_NAAPS_compare_types(date_str, minlat = 65., ax = None, \
        vmin = None, vmax = None, plot_log = False, ptitle = '', \
        circle_bound = True, zoom = True, \
        save = save_flag)

date_str = '20080422'
plot_NAAPS_compare_types(date_str, minlat = 65., ax = None, \
    vmin = None, vmax = None, plot_log = False, ptitle = '', \
    circle_bound = True, zoom = True, \
    save = save_flag)

sys.exit()

NAAPS_noAI   = read_NAAPS_event('20080422', minlat = minlat, dtype = 'no_AI')
NAAPS_withAI = read_NAAPS_event('20080422', minlat = minlat, dtype = 'with_AI')

sys.exit()

# NOTE: Remake gridded NAAPS data to include
#       drysink and wetsink?

NAAPS_data = readgridNAAPS_NCDF(infile=home_dir + \
    '/Research/NAAPS/naaps_grid_smoke_conc_sfc_2005_2020.nc',\
    start_date = 200504, end_date = 202009, calc_month = True, \
    minlat = 65)

plotNAAPS_ClimoTrend_all(NAAPS_data,\
    trend_type = 'standard', minlat=70,save=True)

sys.exit()

#date_str = '20080422'
date_str = '20120615'
#date_str = '20140726'
#date_str = '20140802'
#date_str = '20170816'
#date_str = '20200824'
#var = 'smoke_wetsink'
var = 'smoke_conc_sfc'

if(date_str == '20080422'):
    minlat = 60.
    min_smoke = 20
    max_smoke = 2e5
    vmax = 300
    lat_bounds = [65, 75]
    lon_bounds = [150, 210]
    min_ice = 80.
    vmin2 = 0.4
elif(date_str == '20120615'):
    begin_str = '20120501'
    end_str   = '20120731'
    #begin_str = '20120501'
    #end_str   = '20120731'
    interval = 4 
    minlat = 60.
    min_smoke = 15
    max_smoke = 35
    #max_smoke = 50
    vmax = 50
    lat_bounds = [60, 75]
    lon_bounds = [305, 345]
    min_ice = 0.
    alb_min = 0.4
    ##!#plats =  [\
    ##!#    64.5616,  # smoke
    ##!#    66.6799,  # smoke
    ##!#    62.4603,
    ##!#    69.7040,  # light smoke
    ##!#    70.3616,  # no smoke
    ##!#    73.6032,  # no smoke
    ##!#]
    ##!#plons =  [\
    ##!#    313.3088,  # smoke
    ##!#    312.6112,  # smoke
    ##!#    314.3501, 
    ##!#    314.7316,  # light smoke
    ##!#    327.0151,  # no smoke
    ##!#    316.2099,  # no smoke
    ##!#]
elif(date_str == '20140802'):
    begin_str = '20140701'
    end_str   = '20140831'
    interval = 4 
    minlat = 60.
    min_smoke = 30
    max_smoke = 150
    vmax = 50
    lat_bounds = [58, 82]
    lon_bounds = [280, 355]
    min_ice = 0.
    alb_min = 0.4
if(date_str == '20170816'):
    begin_str = '20170801'
    end_str   = '20170831'
    interval = 4
    minlat = 72.
    #minlat = 78.
    min_smoke = 10
    max_smoke = 300
    vmax = 300
    lat_bounds = [minlat, 83]
    #lat_bounds = [minlat, 85]
    lon_bounds = [290, 340]
    #lon_bounds = [180, 340]
    min_ice = 70.
    alb_min = 0.2
elif(date_str == '20200824'):
    minlat = 75.
    min_smoke = 10
    max_smoke = 150
    vmax = 150
    lat_bounds = [75, 90]
    lon_bounds = [0, 360]
    min_ice = 0.
    vmin2 = 0.2


dt_begin_local1 = datetime(2012,6,15)
dt_end_local1 = datetime(2012,6,19)
NAAPS_data = read_NAAPS_event('20120615', minlat = minlat)
NCEP_data = read_zoom_NCEP_data_single(\
    dt_begin_local1.strftime('%Y%m%d'), \
    dt_end_local1.strftime('%Y%m%d'), \
    min_ice, \
    min_smoke, max_smoke, alb_min, NAAPS_data, minlat, lat_bounds, \
    lon_bounds, satellite = 'All', mask_NAAPS = False, \
    plot_daily_data = True)

sys.exit()

plot_NCEP_region_time_series_combined(date_str, begin_str, end_str, \
    interval, var, \
    minlat = minlat, vmin = None, vmax = vmax, vmin2 = alb_min, vmax2 = None, \
    min_ice = min_ice, min_smoke = min_smoke, max_smoke = max_smoke, plot_log = False, \
    satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
    lon_bounds = lon_bounds, ax = None, plot_daily_data = False, \
    zoom = True, save = True)

sys.exit()



plot_NAAPS_event_CERES_region_time_series_combined(date_str, begin_str, end_str, 
    interval, var, ceres_var = 'lwf_clr', minlat = minlat, vmin = None, \
    vmax = vmax, vmin2 = alb_min, vmax2 = None, min_ice = min_ice, \
    min_smoke = min_smoke, max_smoke = max_smoke, plot_log = False, \
    satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
    lon_bounds = lon_bounds, plot_daily_data = False, \
    zoom = True, save = True)


sys.exit()



plot_MODIS_data_before_after(date_str, save = False)
sys.exit()




plot_NAAPS_event_CERES_region_comp_twovars(date_str, var, ceres_var1 = 'swf_clr', \
    ceres_var2 = 'lwf_clr', \
    minlat = minlat, vmin = None, vmax = vmax, vmin2 = alb_min, vmax2 = None, \
    min_ice = 80., min_smoke = 0, max_smoke = max_smoke, plot_log = True, \
    satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
    lon_bounds = lon_bounds, plot_daily_data = False, \
    zoom = True, save = True)

sys.exit()

#plot_MODIS_OMI_data(date_str, save = True)
#


#plot_NAAPS_event_CERES_region_time_series_combined(date_str, begin_date, end_date, \
#    interval, var, ceres_var = 'alb_clr', \
#    minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
#    min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
#    satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
#    lon_bounds = [-180, 180], plot_daily_data = False, \
#    zoom = True, save = False):

#combined_data, combined_data_nosmoke = \
plot_NAAPS_CERES_flux_diffs(date_str, var, ceres_var = 'alb_clr', \
    minlat = minlat, vmin = None, vmax = vmax, vmin2 = alb_min, vmax2 = None, \
    min_ice = 80., min_smoke = 0, max_smoke = max_smoke, plot_log = True, \
    satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
    lon_bounds = lon_bounds, plot_daily_data = False, \
    in_year = None, markersize = None, zoom = True, save = True)

sys.exit()


pvars = ['alb_clr', 'swf_clr', 'lwf_clr']
for cvar in pvars:
    plot_NAAPS_multi_CERES_region_time_series(date_str, begin_str, end_str, 
        interval, var, ceres_var = cvar, minlat = minlat, vmin = None, \
        vmax = vmax, vmin2 = alb_min, vmax2 = None, min_ice = min_ice, \
        min_smoke = min_smoke, max_smoke = max_smoke, plot_log = False, \
        satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
        lon_bounds = lon_bounds, plot_daily_data = False, \
        zoom = True, save = True)

sys.exit()

# beg: 2020082400
# end: 2020090800
base_date = datetime.strptime('2012070200','%Y%m%d%H')
for ii in range(60):
    new_date = base_date + timedelta(hours = ii * 6)
    plot_NAAPS_figure(new_date.strftime('%Y%m%d%H'), var, minlat = 60., \
        vmax = 5, plot_log = False, zoom = True, save = True)

sys.exit()


plot_NAAPS_event_CERES_region_time_series_allvars(date_str, begin_str, end_str, 
    interval, var, minlat = minlat, vmin = None, \
    vmax = vmax, vmin2 = alb_min, vmax2 = None, min_ice = min_ice, \
    min_smoke = min_smoke, max_smoke = max_smoke, plot_log = False, \
    satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
    lon_bounds = lon_bounds, plot_daily_data = False, \
    zoom = True, save = False)

sys.exit()

second_arrays, second_labels, second_arrays_nosmoke = \
    plot_NAAPS_event_CERES(date_str, var, ceres_var = 'swf_clr', \
#ceres = plot_NAAPS_event_CERES(date_str, var, ceres_var = 'alb_clr', \
    satellite = 'All', minlat = minlat, vmin = None, vmax = vmax, \
    #satellite = 'All', minlat = 60., vmin = None, vmax = 300, \
    vmin2 = alb_min, vmax2 = 300, \
    plot_log = False, min_ice = min_ice, min_smoke = min_smoke, max_smoke = max_smoke, \
    ptitle = '', zoom = True, lat_bounds = lat_bounds, lon_bounds = lon_bounds, \
    plot_daily_data = False, save = True)
sys.exit()


#plot_NAAPS_event_CERES_region_comp(date_str, var, ceres_var = 'alb_clr', \
#    minlat = minlat, vmin = None, vmax = vmax, vmin2 = alb_min, vmax2 = None, \
#    min_ice = 80., min_smoke = 0, max_smoke = max_smoke, plot_log = True, \
#    satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
#    lon_bounds = lon_bounds, plot_daily_data = False, \
#    zoom = True, save = False)

pvars = ['alb_clr', 'swf_clr', 'lwf_clr']
for cvar in pvars:
    #plot_NAAPS_event_CERES_region_comp(date_str, var, ceres_var = cvar, \
    #    minlat = minlat, vmin = None, vmax = vmax, vmin2 = alb_min, vmax2 = None, \
    #    min_ice = 80., min_smoke = 0, max_smoke = max_smoke, plot_log = True, \
    #    satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
    #    lon_bounds = lon_bounds, plot_daily_data = False, \
    #    zoom = True, save = True)
    second_arrays, second_labels, second_arrays_nosmoke = \
            plot_NAAPS_event_CERES_region_time_series(date_str, begin_str, end_str, 
            interval, var, ceres_var = cvar, minlat = minlat, vmin = None, \
            vmax = vmax, vmin2 = alb_min, vmax2 = None, min_ice = min_ice, \
            min_smoke = min_smoke, max_smoke = max_smoke, plot_log = False, \
            satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
            lon_bounds = lon_bounds, plot_daily_data = False, \
            zoom = True, save = True)
sys.exit()

sys.exit()




#plot_NAAPS_multi_CERES_region_comp(date_str, var, ceres_var = 'lwf_clr', \

plot_NAAPS_event(date_str, var, minlat = 60., vmin = None, vmax = 200,\
    plot_log = False, ptitle = '', zoom = True, save = False)
sys.exit()


plot_NAAPS_event_CERES_points_time_series(date_str, begin_str, end_str, \
    interval, var, plats, plons, ceres_var = 'alb_clr', \
    minlat = minlat, vmin = None, vmax = vmax, vmin2 = alb_min, vmax2 = None, \
    min_ice = min_ice, min_smoke = min_smoke, max_smoke = max_smoke, \
    plot_log = False, \
    satellite = 'All', ptitle = '', lat_bounds = lat_bounds, \
    lon_bounds = lon_bounds, plot_daily_data = False, \
    zoom = True, save = True)

sys.exit()

NAAPS_data = readgridNAAPS_NCDF(infile=home_dir + \
    '/Research/NAAPS/naaps_grid_smoke_conc_sfc_2005_2020.nc',\
    start_date = 200504, end_date = 202009, calc_month = True, \
    minlat = 65)

plotNAAPS_ClimoTrend_all(NAAPS_data,\
    trend_type = 'standard', minlat=65,save=True)

sys.exit()

#plotNAAPS_MonthTrend(NAAPS_data,month_idx=4,save=False,\
#    trend_type='standard',minlat=65.,return_trend=False, colorbar = True, \
#    title = '', label = '', colorbar_label_size = 14, pax = None, \
#    show_pval = True, uncert_ax = None)

sys.exit()

# NOTE: to regenerate the entire NAAPS monthly climatology dataset
##!#NAAPS_data = calc_NAAPS_all_avgs('200504', '202009', minlat = 65., \
##!#    mask_zero = False)
##!#sys.exit()





# Bugged
#2013 466 501
#2014 488 500
#2015 452 460
#2016 392 477
#2017 249 505
#2018 408 488
#2019 452 452
#2020 464 464
#2021 465 484

# Fixed
#2013 466 489
#2014 488 485
#2015 452 451
#2016 392 476
#2017 249 391
#2018 408 452
#2019 452 243
#2020 464 356
#2021 465 391


for ii in range(second_arrays.shape[0]):
    u_stat, p_val = mannwhitneyu(second_arrays[ii,0], second_arrays[ii,1],\
        alternative = 'greater')
    print(ii, p_val)

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#ax.hist(second_arrays[4,1], label = 'After')
#ax.hist(second_arrays[4,0], label = 'Before')
#ax.legend()
plt.show()


